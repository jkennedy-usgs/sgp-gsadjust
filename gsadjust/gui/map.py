from PyQt5 import QtCore, QtWidgets
from functools import partial
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.dates import date2num
from matplotlib.figure import Figure

from .section import Section
import numpy as np
import datetime as dt
import logging
from .messages import MessageBox
from ..data.nwis import get_gwl_change, plot_hydrograph_with_gravity


class GravityChangeMap(QtWidgets.QDialog):
    """
    Map window for showing gravity-change results

    Parameters
    ----------
    table
    header
    coords
    full_table
    full_header

    """

    station_label = None
    try:
        import cartopy.crs as ccrs
        import cartopy.io.img_tiles as cimgt
        from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
        from cartopy.mpl.ticker import LongitudeLocator, LatitudeLocator
    except ImportError as e:
        logging.info(e)

    def __init__(self, table, header, coords, full_table, full_header, parent=None):
        super(GravityChangeMap, self).__init__(parent)
        self.plot_buttons = []
        # table: pre-calculated survey-to-survey change
        self.table = [list(i) for i in zip(*table)]
        self.header = header
        self.coords = coords
        # full_table: for calculating change between any 2 surveys
        self.full_table = full_table  # [list(i) for i in zip(*full_table)]
        self.full_header = full_header
        # instance variables
        self.figure = Figure(figsize=(10, 8))
        self.n_surveys = (len(self.table) - 1) / 2
        self.surveys = self.get_survey_dates(header)
        self.axlim = self.get_axis_lims_from_data(coords)
        self.default_axlim = self.axlim
        self.g_cb = None  # matplotlib colorbar
        self.wl_cb = None
        self.points = None  # matplotlib PathCollection, returned from scatter
        self.wl_points = None  #
        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvasLayout = QtWidgets.QVBoxLayout()
        self.canvas_widget = QtWidgets.QWidget()
        # self.canvas_widget.resize(400,400)
        self.canvas = FigureCanvas(self.figure)
        self.canvasLayout.addWidget(self.canvas, stretch=1)
        self.canvas_widget.setLayout(self.canvasLayout)
        # self.canvas.setParent(self.canvas_widget)
        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.home = self.home_button

        # Choose time interval
        self.time_widget = QtWidgets.QWidget()
        self.btnIncremental = QtWidgets.QRadioButton("Incremental", self)
        self.btnIncremental.setChecked(True)
        self.btnReference = QtWidgets.QRadioButton("Change relative to", self)
        self.drpReference = QtWidgets.QComboBox()
        self.btnTrend = QtWidgets.QRadioButton("Trend", self)
        self.btnIncremental.toggled.connect(self.plot)
        self.btnReference.toggled.connect(self.plot)
        self.drpReference.currentIndexChanged.connect(self.plot)
        self.btnTrend.toggled.connect(self.plot)
        bbox1 = QtWidgets.QHBoxLayout()
        bbox1.setContentsMargins(0, 0, 0, 0)
        bbox1.addWidget(self.btnIncremental)
        bbox1.addSpacing(10)
        bbox1.addWidget(self.btnReference)
        bbox1.addWidget(self.drpReference)
        bbox1.addSpacing(10)
        bbox1.addWidget(self.btnTrend)
        bbox1.addStretch(1)
        self.time_widget.setLayout(bbox1)

        self.basemap_widget = QtWidgets.QWidget()
        self.cbBasemap = QtWidgets.QCheckBox("Show Basemap", self)
        self.cbBasemap.setChecked(False)
        self.cbBasemap.stateChanged.connect(self.plot)
        self.drpBasemap = QtWidgets.QComboBox()
        self.drpBasemap.currentIndexChanged.connect(self.plot)
        self.sliderResolution = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.sliderResolution.setMinimum(5)
        self.sliderResolution.setMaximum(17)
        self.sliderResolution.setValue(12)
        btnRefresh = QtWidgets.QPushButton("Refresh")
        btnRefresh.clicked.connect(self.resolution_slider_changed)
        bbox2 = QtWidgets.QHBoxLayout()
        bbox2.setContentsMargins(0, 0, 0, 0)
        bbox2.addWidget(self.cbBasemap)
        bbox2.addSpacing(10)
        bbox2.addWidget(self.drpBasemap)
        bbox2.addSpacing(40)
        bbox2.addStretch(1)
        bbox2.addWidget(QtWidgets.QLabel("Zoom"))
        bbox2.addWidget(self.sliderResolution)
        bbox2.addSpacing(10)
        bbox2.addWidget(btnRefresh)

        self.basemap_widget.setLayout(bbox2)

        self.units_widget = QtWidgets.QWidget()
        self.cbUnits = QtWidgets.QCheckBox("Show change in meters of water", self)
        self.cbUnits.setChecked(False)
        self.cbUnits.stateChanged.connect(self.plot)
        self.sliderColorRange = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.sliderColorRange.setMinimum(1)
        self.sliderColorRange.setMaximum(15)
        self.sliderColorRange.setValue(100)
        self.sliderColorRange.resize(100, 10)
        self.sliderColorRange.setTickInterval(10)
        self.sliderColorRange.valueChanged.connect(self.update_plot)
        bbox = QtWidgets.QHBoxLayout()
        bbox.setContentsMargins(0, 0, 0, 0)
        bbox.addWidget(self.cbUnits)
        bbox.addSpacing(40)
        bbox.addStretch(1)
        bbox.addWidget(QtWidgets.QLabel("Color range"))
        bbox.addWidget(self.sliderColorRange)
        bbox.addSpacing(80)
        self.units_widget.setLayout(bbox)

        self.waterlevels_widget = QtWidgets.QWidget()
        self.cbWaterlevels = QtWidgets.QCheckBox("Show NWIS water levels")
        self.cbWaterlevels.setChecked(False)
        self.cbWaterlevels.stateChanged.connect(self.plot)
        self.gwl_sliderColorRange = QtWidgets.QSlider(QtCore.Qt.Horizontal)

        self.gwl_sliderColorRange.setMinimum(1)
        self.gwl_sliderColorRange.setMaximum(10)
        self.gwl_sliderColorRange.setValue(2)
        self.gwl_sliderColorRange.resize(100, 10)
        self.gwl_sliderColorRange.valueChanged.connect(self.update_plot)
        self.gwl_datewindow = QtWidgets.QSpinBox()

        btnRefresh2 = QtWidgets.QPushButton("Refresh")
        btnRefresh2.clicked.connect(self.plot)
        bbox3 = QtWidgets.QHBoxLayout()
        bbox3.setContentsMargins(0, 0, 0, 0)
        bbox3.addWidget(self.cbWaterlevels)
        bbox3.addSpacing(40)
        bbox3.addWidget(QtWidgets.QLabel("Date window (days)"))
        bbox3.addWidget(self.gwl_datewindow)
        bbox3.addSpacing(40)
        bbox3.addStretch(1)
        bbox3.addWidget(QtWidgets.QLabel("Color range"))
        bbox3.addWidget(self.gwl_sliderColorRange)
        bbox3.addWidget(btnRefresh2)
        self.waterlevels_widget.setLayout(bbox3)

        # Right-side controls (second tuple item is the default)
        self.ctrl_items = {
            "Style": ("Powerpoint", "Print"),
            "Units": ("Miles", "Km"),
            "Tick labels": ("On", "Off"),
            "Scale bar": ("On", "Off"),
            "Colorbar": ("Bottom", "Right"),
        }

        controls_right_widget = QtWidgets.QWidget()
        controls_right_layout = QtWidgets.QVBoxLayout()
        for i in self.ctrl_items.keys():
            item_widget = QtWidgets.QWidget()
            item_layout = QtWidgets.QHBoxLayout()
            item_layout.setContentsMargins(0, 0, 0, 0)
            item_layout.addWidget(QtWidgets.QLabel(i))
            item_layout.addStretch(1)
            r1 = QtWidgets.QRadioButton(self.ctrl_items[i][0])
            r2 = QtWidgets.QRadioButton(self.ctrl_items[i][1])
            r1.setFixedWidth(80)
            r2.setFixedWidth(60)
            r2.setChecked(True)
            r1.clicked.connect(self.plot)
            r2.clicked.connect(self.plot)
            setattr(self, i, (r1, r2))
            item_layout.addWidget(r1)
            item_layout.addWidget(r2)
            item_widget.setLayout(item_layout)
            item_widget.setFixedHeight(20)
            controls_right_layout.addWidget(item_widget)
        controls_right_layout.addStretch(1)
        controls_right_widget.setLayout(controls_right_layout)

        # Date slider
        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider.setTickPosition(QtWidgets.QSlider.TicksBelow)

        # Groundwater-level table
        nwis_section = Section("Groundwater Levels", 100, self)

        nwis_layout = QtWidgets.QVBoxLayout()
        self.gwl_table = QtWidgets.QTableWidget()
        self.gwl_table.setColumnCount(7)
        self.gwl_table.itemClicked.connect(self.gwl_table_checked)
        nwis_layout.addWidget(self.gwl_table)
        nwis_section.setContentLayout(nwis_layout)
        self.slider_label = QtWidgets.QLabel(self.header[1])
        # set the layout
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)

        controls_widget = QtWidgets.QWidget()
        controls_layout = QtWidgets.QHBoxLayout()

        controls_left_widget = QtWidgets.QWidget()
        controls_left_layout = QtWidgets.QVBoxLayout()
        controls_left_layout.addWidget(self.time_widget)
        controls_left_layout.addWidget(self.basemap_widget)
        controls_left_layout.addWidget(self.units_widget)
        controls_left_layout.addWidget(self.waterlevels_widget)
        controls_left_widget.setLayout(controls_left_layout)

        controls_layout.addWidget(controls_left_widget)
        divider = QtWidgets.QFrame()
        divider.setFrameShape(QtWidgets.QFrame.VLine)
        controls_layout.addWidget(divider)
        controls_layout.addWidget(controls_right_widget)
        controls_widget.setLayout(controls_layout)
        #
        # controls_left_widget.resize(200, 20)
        # controls_widget.resize(600,20)

        layout.addWidget(controls_widget, 0)

        layout.addWidget(self.canvas_widget, 1)
        layout.addWidget(self.slider_label)

        self.setWindowTitle("Gravity Change Map")
        layout.addWidget(self.slider)
        layout.addWidget(nwis_section)
        self.setLayout(layout)

        for s in self.surveys:
            self.drpReference.addItem(s)

        self.maps = {
            0: ("ESRI", "NatGeo_World_Map"),
            1: ("ESRI", "USA_Topo_Maps"),
            2: ("ESRI", "World_Imagery"),
            3: ("ESRI", "World_Shaded_Relief"),
            4: ("ESRI", "World_Street_Map"),
            5: ("ESRI", "World_Topo_Map"),
            6: ("USGS", "USGSHydroCached"),
            7: ("USGS", "USGSImageryOnly"),
            8: ("USGS", "USGSImageryTopo"),
            9: ("USGS", "USGSShadedReliefOnly"),
            10: ("USGS", "USGSTopo"),
        }
        for m in self.maps.values():
            self.drpBasemap.addItem(m[1])
        self.plot()

    def home_button(self):
        # This isn't being called, but Home button seems to work okay? Need to figure
        # out which instances the home button isn't working.
        self.axlim = self.default_axlim
        self.plot()

    def resolution_slider_changed(self):
        if self.cbBasemap.isChecked():
            self.plot()

    def get_survey_dates(self, header):
        dates = []
        first = True
        for title in header:
            if first:
                first = False
                continue
            d = title.split("_to_")
            dates.append(d[0])
            dates.append(d[1])
        return sorted(list(set(dates)))

    def update_plot(self):

        try:
            self.slider.setEnabled(False)
            self.plot_g()
            if self.cbWaterlevels.isChecked():
                self.plot_wls()
            # else:
            #     try:
            #         self.wl_cb.remove()
            #     except (KeyError, AttributeError) as e:
            #         pass
            self.slider_label.setText(self.get_name())
            if getattr(self, "Style")[0].isChecked():  # Powerpoint
                self.ax.set_title(self.get_name(), fontsize=16, fontweight="bold")
            else:
                self.ax.set_title(self.get_name(), fontsize=12)
            self.canvas.draw()
            self.slider.setEnabled(True)
        except TypeError as e:
            # Probably a missing coordinate
            j = 1
            return

    def plot_g(self):
        try:
            # For some reason I can't figure out, update_plot is called multiple
            # times when changing the slider. The .remove() calls cause exceptions
            # on the subsequent calls because they've already been removed.
            self.g_cb.remove()
            self.points.remove()
        except (KeyError, ValueError) as e:
            pass
        except AttributeError as e:
            pass
        try:
            x, y, d, names = self.get_g_data()
        except KeyError as e:
            MessageBox.warning(
                "Map view gravity change",
                f"Could not plot station {e.args[0]}. Check that station names in the "
                "coordinates table match the names in the survey.",
            )
            return
        clim = self.get_color_lims()

        self.points = self.ax.scatter(
            x,
            y,
            c=d,
            s=200,
            vmin=clim[0],
            vmax=clim[1],
            cmap="RdBu",
            picker=5,
            zorder=10,
        )
        if getattr(self, "Colorbar")[0].isChecked():  # Bottom
            self.g_cb = self.figure.colorbar(
                self.points, orientation="horizontal", ax=self.ax
            )
        else:
            self.g_cb = self.figure.colorbar(self.points, ax=self.ax)
        # self.g_cb = self.figure.colorbar(self.points, ax=self.ax)
        self.set_g_cb_label()
        self.points.names = names
        # self.canvas.draw()

    def plot_wls(self):
        try:
            # For some reason I can't figure out, update_plot is called multiple
            # times when changing the slider. The .remove() calls cause exceptions
            # on the subsequent calls because they've already been removed.
            self.wl_points.remove()
            self.wl_cb.remove()
        except (KeyError, ValueError) as e:
            self.wl_cb = None
            self.wl_points = None
            pass
        except AttributeError as e:
            self.wl_cb = None
            self.wl_points = None
            pass
        clim = self.get_gwl_color_lims()
        gwl_x, gwl_y, gwl_d, names = self.get_checked_gwl_data()
        self.wl_points = self.ax.scatter(
            gwl_x,
            gwl_y,
            c=gwl_d,
            s=50,
            marker="s",
            edgecolors="k",
            vmin=clim[0],
            vmax=clim[1],
            cmap="BrBG",
            picker=5,
            zorder=10,
        )
        if self.wl_cb is None:
            if getattr(self, "Colorbar")[0].isChecked():  # Bottom
                self.wl_cb = self.figure.colorbar(
                    self.wl_points, orientation="horizontal"
                )
            else:
                self.wl_cb = self.figure.colorbar(self.wl_points)
        self.wl_cb.mappable.set_clim(vmin=clim[0], vmax=clim[1])
        self.set_wl_cb_label()
        self.wl_points.names = names

    def get_datenums(self, full_header):
        # Get dates as datetimes for trend analysis
        obs_idxs, obs_dates = [], []
        for item in full_header:
            if item[-2:] == "_g":
                try:
                    obs_dates.append(
                        date2num(dt.datetime.strptime(item[:-2], "%Y-%m-%d"))
                    )
                    obs_idxs.append(full_header.index(item))
                except TypeError as e:
                    j = 1
                    pass
        return obs_dates, obs_idxs

    def get_checked_gwl_data(self):
        x, y, d, names = [], [], [], []
        for i in range(self.gwl_table.rowCount()):
            item = self.gwl_table.item(i, 0)
            if item.checkState() == QtCore.Qt.Checked:
                x.append(float(self.gwl_table.item(i, 2).text()))
                y.append(float(self.gwl_table.item(i, 1).text()))
                d.append(float(self.gwl_table.item(i, 3).text()))
                names.append(self.gwl_table.item(i, 0).text())
        return x, y, d, names

    def get_g_data(self):
        x, y, value, name = [], [], [], []

        if self.btnIncremental.isChecked():
            data = self.table[self.slider.value()]
            stations = self.table[0]
            for sta_idx, sta in enumerate(stations):
                try:
                    datum = float(data[sta_idx])
                except ValueError as e:
                    continue
                try:
                    x.append(self.coords[sta][0])
                    y.append(self.coords[sta][1])
                    if self.cbUnits.isChecked():
                        datum /= 41.9
                    value.append(datum)
                    name.append(sta)
                except KeyError as e:
                    MessageBox.critical(
                        "Error", f"Coordinates for station {sta} not found"
                    )
                    return None

        elif self.btnReference.isChecked():
            ref_survey = (
                self.drpReference.currentData(role=QtCore.Qt.DisplayRole) + "_g"
            )
            ref_col_idx = self.full_header.index(ref_survey)
            current_survey = self.surveys[self.slider.value() - 1] + "_g"
            current_col_idx = self.full_header.index(current_survey)
            for r in self.full_table:
                sta = r[0]
                try:
                    ref_g = float(r[ref_col_idx])
                    surv_g = float(r[current_col_idx])
                except ValueError as e:
                    continue
                x.append(self.coords[sta][0])
                y.append(self.coords[sta][1])

                if current_survey < ref_survey:
                    datum = ref_g - surv_g
                else:
                    datum = surv_g - ref_g
                if self.cbUnits.isChecked():
                    value.append(datum / 41.9)
                else:
                    value.append(datum)
                name.append(sta)

        elif self.btnTrend.isChecked():
            obs_dates, obs_idxs = self.get_datenums(self.full_header)
            for r in self.full_table:
                sta = r[0]
                X = []
                Y = [float(r[idx]) for idx in obs_idxs if float(r[idx]) > -998]
                for idx, obs_idx in enumerate(obs_idxs):
                    if float(r[obs_idx]) > -998:
                        X.append(obs_dates[idx])
                if len(X) < 1:
                    continue
                z = np.polyfit(X, Y, 1)
                x.append(self.coords[sta][0])
                y.append(self.coords[sta][1])
                uGal_per_year = z[0] * 365.25
                if self.cbUnits.isChecked():
                    value.append(uGal_per_year / 41.9)
                else:
                    value.append(uGal_per_year)
                name.append(sta)

        return x, y, value, name

    def get_name(self):
        if self.btnTrend.isChecked():
            title = "{} to {}".format(self.surveys[0], self.surveys[-1])
            return title
        elif self.btnIncremental.isChecked():
            name = self.header[self.slider.value()]
            return name.replace("_", " ")
        elif not self.btnIncremental.isChecked():
            ref_survey = self.drpReference.currentData(role=QtCore.Qt.DisplayRole)
            current_survey = self.surveys[self.slider.value() - 1]
            if dt.datetime.strptime(current_survey, "%Y-%m-%d") < dt.datetime.strptime(
                ref_survey, "%Y-%m-%d"
            ):
                return "{} to {}".format(current_survey, ref_survey)
            else:
                return "{} to {}".format(ref_survey, current_survey)

    def plot(self):
        if self.btnTrend.isChecked():
            self.slider.setEnabled(False)
        if self.btnIncremental.isChecked():
            self.slider.setEnabled(True)
            self.slider.setRange(1, self.n_surveys)
        elif not self.btnIncremental.isChecked():
            self.slider.setEnabled(True)
            self.slider.setRange(1, self.n_surveys + 1)
        self.figure.clf()
        self.ax = self.figure.add_subplot(
            1,
            1,
            1,
            projection=self.ccrs.PlateCarree(),
        )
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        self.plot_g()

        if self.cbWaterlevels.isChecked():
            dates = self.get_name().split(" to ")
            try:
                self.populate_gwl_table(dates)
                self.plot_wls()
            except IndexError as e:
                pass
        if getattr(self, "Style")[0].isChecked():  # Powerpoint
            self.ax.set_title(self.get_name(), fontsize=16, fontweight="bold")
        else:
            self.ax.set_title(self.get_name(), fontsize=12)
        self.ax.set_extent(self.axlim)

        if self.cbBasemap.isChecked():
            self.show_background(self.sliderResolution.value())

        if getattr(self, "Tick labels")[0].isChecked():
            xticks = self.LongitudeLocator(nbins=2)._raw_ticks(
                self.axlim[0], self.axlim[1]
            )
            yticks = self.LatitudeLocator(nbins=3)._raw_ticks(
                self.axlim[2], self.axlim[3]
            )
            self.ax.set_xticks(xticks, crs=self.ccrs.PlateCarree())
            self.ax.set_yticks(yticks, crs=self.ccrs.PlateCarree())
            self.ax.xaxis.set_major_formatter(
                self.LongitudeFormatter(dms=True, auto_hide=False)
            )
            self.ax.yaxis.set_major_formatter(
                self.LatitudeFormatter(dms=True, auto_hide=False)
            )
            self.ax.tick_params(direction="in")
        if getattr(self, "Scale bar")[0].isChecked():
            self.scale_bar(self.ax)
        self.ax.callbacks.connect("xlim_changed", self.on_lims_change)
        self.ax.callbacks.connect("ylim_changed", self.on_lims_change)
        self.slider_label.setText(self.get_name())

        # refresh canvas
        self.figure.canvas.mpl_connect("pick_event", self.show_point_label)
        self.canvas.draw()

        self.slider.valueChanged.connect(self.update_plot)
        QtWidgets.QApplication.restoreOverrideCursor()

    def show_hydrograph(self, station):
        dates = self.get_name().split(" to ")
        plot_hydrograph_with_gravity(station, dates)

    def populate_gwl_table(self, dates):
        win = int(self.gwl_datewindow.text())
        if getattr(self, "Units")[0].isChecked():  # miles
            convert_to_meters = False
        else:
            convert_to_meters = True

        sites, nwis_dg, lat, lon, date1, date2 = get_gwl_change(
            dates, self.axlim, buffer=win, convert=convert_to_meters
        )
        self.gwl_table.setRowCount(len(sites))
        for idx, station in enumerate(sites):
            item = QtWidgets.QTableWidgetItem(station)
            item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
            item.setCheckState(QtCore.Qt.Checked)
            btn = QtWidgets.QPushButton("Hydrograph")
            btn.clicked.connect(partial(self.show_hydrograph, station))
            self.plot_buttons.append(btn)
            self.gwl_table.setItem(idx, 0, item)
            self.gwl_table.setItem(
                idx, 1, QtWidgets.QTableWidgetItem(f"{lat[station]:.4f}")
            )
            self.gwl_table.setItem(
                idx, 2, QtWidgets.QTableWidgetItem(f"{lon[station]:.4f}")
            )
            self.gwl_table.setItem(
                idx, 3, QtWidgets.QTableWidgetItem(f"{nwis_dg[station]:.2f}")
            )
            self.gwl_table.setItem(
                idx, 4, QtWidgets.QTableWidgetItem(f"{date1[station]}")
            )
            self.gwl_table.setItem(
                idx, 5, QtWidgets.QTableWidgetItem(f"{date2[station]}")
            )
            self.gwl_table.setCellWidget(idx, 6, self.plot_buttons[-1])
            if getattr(self, "Units")[0].isChecked():  # miles
                self.gwl_table.setHorizontalHeaderLabels(
                    [
                        "Station",
                        "Latitude",
                        "Longitude",
                        "GWL change (ft)",
                        "From",
                        "To",
                        "",
                    ]
                )
            else:
                self.gwl_table.setHorizontalHeaderLabels(
                    [
                        "Station",
                        "Latitude",
                        "Longitude",
                        "GWL change (m)",
                        "From",
                        "To",
                        "",
                    ]
                )
            self.gwl_table.setColumnWidth(0, 120)

    def gwl_table_checked(self, item):
        return
        # # self.update_plot()
        # if item.checkState() == QtCore.Qt.Checked:
        #     self.update_plot()
        # else:
        #     print('"%s" Clicked' % item.text())

    def set_g_cb_label(self):
        if getattr(self, "Style")[0].isChecked():  # Powerpoint
            fs = 16
        else:
            fs = 12
        if self.btnTrend.isChecked():
            if not self.cbUnits.isChecked():
                self.g_cb.set_label("Gravity trend, in µGal/yr", fontsize=fs)
            else:
                self.g_cb.set_label("Gravity trend, in meters of water/yr", fontsize=fs)
        else:
            if not self.cbUnits.isChecked():
                self.g_cb.set_label("Gravity change, in µGal", fontsize=fs)
            elif self.cbUnits.isChecked():
                self.g_cb.set_label(
                    "Aquifer-storage change,\n in meters of water", fontsize=fs
                )

    def set_wl_cb_label(self):
        if getattr(self, "Style")[0].isChecked():  # Powerpoint
            fs = 16
        else:
            fs = 12
        if getattr(self, "Units")[0].isChecked():  # P
            self.wl_cb.set_label("Groundwater-level change, in feet", fontsize=fs)
        else:
            self.wl_cb.set_label("Groundwater-level change, in meters", fontsize=fs)

    def on_lims_change(self, axes):
        self.axlim = self.ax.get_extent()
        return

    def show_background(self, zoom):
        server, m = self.maps[self.drpBasemap.currentIndex()]
        if server == "ESRI":
            self.stamen_terrain = self.cimgt.GoogleTiles(
                url="https://server.arcgisonline.com/ArcGIS/rest/services/"
                + m
                + "/MapServer/tile/{z}/{y}/{x}.jpg"
            )
        # USGS server omits the .jpg
        else:
            self.stamen_terrain = self.cimgt.GoogleTiles(
                url="https://basemap.nationalmap.gov/arcgis/rest/services/"
                + m
                + "/MapServer/tile/{z}/{y}/{x}"
            )

        self.ax.add_image(
            self.stamen_terrain, zoom, interpolation="spline36", regrid_shape=2000
        )

    def show_point_label(self, event):
        """
        Shows the station name in the upper left of the drift plot when a
        line is clicked.

        Parameters
        ----------
        event : Matplotlib event
        axes : Current axes
            Differs for none|netadj|roman vs continuous

        """

        thispoint = event.artist
        station_label = thispoint.names[event.ind[0]]
        old_station_label = ""
        if self.station_label is not None:
            old_station_label = self.station_label.get_text()
            self.station_label.set_text("")
        if old_station_label != station_label or self.station_label is None:
            self.station_label = self.ax.text(
                0.05,
                0.95,
                station_label,
                horizontalalignment="left",
                verticalalignment="center",
                transform=self.ax.transAxes,
            )
        self.canvas.draw()

        if len(station_label) == 15:
            row = self.find_corresponding_table_row(station_label)
            self.gwl_table.selectRow(row)

    def find_corresponding_table_row(self, station_label):
        for r in range(self.gwl_table.rowCount()):
            if self.gwl_table.item(r, 0).text() == station_label:
                return r

    def get_color_lims(self):
        slider_val = self.sliderColorRange.value() * 10
        clim = (slider_val * -1, slider_val)
        if self.cbUnits.isChecked():
            clim = (clim[0] / 41.9, clim[1] / 41.9)
        return clim

    def get_gwl_color_lims(self):
        slider_val = self.gwl_sliderColorRange.value()
        clim = (slider_val * -0.5, slider_val * 0.5)
        return clim

    def get_axis_lims_from_data(self, coords):
        ratio = 1.2  # width to height
        margin = 0.25

        X, Y, z = zip(*coords.values())

        # Remove values far from the mean (e.g., 0,0)
        n = 2
        X = [x for x in X if abs(x - np.mean(X)) < np.std(X) * n]
        Y = [x for x in Y if abs(x - np.mean(Y)) < np.std(Y) * n]

        xrange = abs(max(X) - min(X))
        yrange = abs(max(Y) - min(Y))

        if xrange > yrange:
            yrange = xrange / ratio
        else:
            xrange = yrange / ratio

        xmin = min(X) - xrange * margin
        xmax = max(X) + xrange * margin
        ymin = min(Y) - yrange * margin
        ymax = max(Y) + yrange * margin

        return xmin, xmax, ymin, ymax

    def scale_bar(
        self,
        ax,
        length=None,
        location=(0.5, 0.07),
        linewidth=1,
    ):
        """
        ax is the axes to draw the scalebar on.
        length is the length of the scalebar in km.
        location is center of the scalebar in axis coordinates.
        (ie. 0.5 is the middle of the plot)
        linewidth is the thickness of the scalebar.
        """
        # Get the limits of the axis in lat long
        llx0, llx1, lly0, lly1 = ax.get_extent(self.ccrs.PlateCarree())
        # Make tmc horizontally centred on the middle of the map,
        # vertically at scale bar location
        sbllx = (llx1 + llx0) / 2
        sblly = lly0 + (lly1 - lly0) * location[1]
        tmc = self.ccrs.TransverseMercator(sbllx, sblly)
        # Get the extent of the plotted area in coordinates in metres
        x0, x1, y0, y1 = ax.get_extent(tmc)
        # Turn the specified scalebar location into coordinates in metres
        sbx = x0 + (x1 - x0) * location[0]
        sby = y0 + (y1 - y0) * location[1]

        # The figure is in degree units, but the scale bar is in meter units, so I
        # couldn't figure out a good way to set a constant length for the scale bar
        # ticks (If specified as a constant, their length would vary depending on the
        # map scale. Here they are set to 1/50 the x-axis length.
        sb_ticklength = (x1 - x0) / 50
        # Calculate a scale bar length if none has been given
        # (Theres probably a more pythonic way of rounding the number but this works)
        if not length:
            length = max([1, (x1 - x0) / 5000])  # in km
            ndim = int(np.floor(np.log10(length)))  # number of digits in number
            length = round(length, -ndim) * 2  # round to 1sf

            # Returns numbers starting with the list
            def scale_number(x):
                if str(x)[0] in ["1", "2", "5"]:
                    return int(x)
                else:
                    return scale_number(x - 10**ndim)

            length = scale_number(length)
        bar_xs_km = [sbx - length * 500, sbx + length * 500]

        if getattr(self, "Style")[0].isChecked():  # Presentation style
            if getattr(self, "Units")[0].isChecked():  # Feet
                length_mi = np.ceil(length * 0.6213)
                bar_xs_mi = [
                    (sbx - length_mi * 500) * 1.6093,
                    (sbx + length_mi * 500) * 1.6093,
                ]
                ax.plot(bar_xs_mi, [sby - 200, sby - 200], transform=tmc, color="k")
                ax.text(
                    sbx,
                    sby - 200,
                    str(length) + " mi",
                    transform=tmc,
                    horizontalalignment="center",
                    verticalalignment="bottom",
                )
            else:  # meters
                ax.plot(bar_xs_km, [sby - 200, sby - 200], transform=tmc, color="k")
                ax.text(
                    sbx,
                    sby - 200,
                    str(length) + " km",
                    transform=tmc,
                    horizontalalignment="center",
                    verticalalignment="bottom",
                )
        else:  # Print style
            # If in feet, the miles scale bar is longer and goes on top. if in m, the km
            # scale bar is longer and goes on top.
            if getattr(self, "Units")[0].isChecked():  # Feet
                length_mi = np.ceil(length * 0.6213)
                bar_xs_primary = [
                    (sbx - length_mi * 500) * 1.6093,
                    (sbx + length_mi * 500) * 1.6093,
                ]
                bar_xs_secondary = bar_xs_km
                # Generate the x coordinate for the ends of the scalebar
                x_locs_btm = [
                    bar_xs_secondary[0],
                    bar_xs_secondary[0]
                    + (bar_xs_secondary[1] - bar_xs_secondary[0]) * 0.25,
                    bar_xs_secondary[0]
                    + (bar_xs_secondary[1] - bar_xs_secondary[0]) * 0.5,
                    bar_xs_secondary[1],
                ]
                x_labels_btm = [0, length * 0.25, length * 0.5, length]
                x_labels_top = [0, length_mi * 0.25, length_mi * 0.5, length_mi]
                # Offset to match the left end with the mile scale bar
                x_locs_btm = [
                    x - (bar_xs_secondary[0] - bar_xs_primary[0]) for x in x_locs_btm
                ]
                x_locs_top = [
                    bar_xs_primary[0],
                    bar_xs_primary[0] + (bar_xs_primary[1] - bar_xs_primary[0]) * 0.25,
                    bar_xs_primary[0] + (bar_xs_primary[1] - bar_xs_primary[0]) * 0.5,
                    bar_xs_primary[1],
                ]
            else:  # Meters; km scale bar on top
                length_mi = np.floor(length * 0.6213)
                bar_xs_primary = [(sbx - length * 500), (sbx + length * 500)]
                bar_xs_secondary = [
                    (sbx - length_mi * 500) * 1.6093,
                    (sbx + length_mi * 500) * 1.6093,
                ]
                x_locs_btm = [
                    bar_xs_secondary[0],
                    bar_xs_secondary[0]
                    + (bar_xs_secondary[1] - bar_xs_secondary[0]) * 0.25,
                    bar_xs_secondary[0]
                    + (bar_xs_secondary[1] - bar_xs_secondary[0]) * 0.5,
                    bar_xs_secondary[1],
                ]
                x_labels_top = [0, length * 0.25, length * 0.5, length]
                x_labels_btm = [0, length_mi * 0.25, length_mi * 0.5, length_mi]
                # Offset to match the left end with the mile scale bar
                x_locs_btm = [
                    x - (bar_xs_secondary[0] - bar_xs_primary[0]) for x in x_locs_btm
                ]
                x_locs_top = [
                    bar_xs_primary[0],
                    bar_xs_primary[0] + (bar_xs_primary[1] - bar_xs_primary[0]) * 0.25,
                    bar_xs_primary[0] + (bar_xs_primary[1] - bar_xs_primary[0]) * 0.5,
                    bar_xs_primary[1],
                ]

            y_locs_down = [sby - sb_ticklength, sby]
            y_locs_up = [sby + sb_ticklength, sby]
            for xl in x_locs_btm:
                ax.plot(
                    [xl, xl], y_locs_down, transform=tmc, color="k", linewidth=linewidth
                )
            for xl in x_locs_top:
                ax.plot(
                    [xl, xl], y_locs_up, transform=tmc, color="k", linewidth=linewidth
                )
            ax.plot(
                bar_xs_primary,
                [sby, sby],
                transform=tmc,
                color="k",
                linewidth=linewidth,
            )

            # Add labels
            text_offset = (
                sb_ticklength + sb_ticklength * 0.1
            )  # Needs to be in meters, relative to the tick length.
            units1 = ""
            units2 = ""
            for idx, x in enumerate(x_labels_top):
                if idx == len(x_labels_top) - 1:
                    if getattr(self, "Units")[0].isChecked():  # Feet
                        units1 = " MI"
                        units2 = " KM"
                    else:
                        units1 = " KM"
                        units2 = " MI"
                if x % 1 < 0.01:
                    x = int(x)
                if x_labels_btm[idx] % 1 < 0.01:
                    x_labels_btm[idx] = int(x_labels_btm[idx])
                ax.text(
                    x_locs_top[idx],
                    sby + text_offset,
                    str(x) + units1,
                    transform=tmc,
                    horizontalalignment="center",
                )
                ax.text(
                    x_locs_btm[idx],
                    sby - text_offset,
                    str(x_labels_btm[idx]) + units2,
                    transform=tmc,
                    verticalalignment="top",
                    horizontalalignment="center",
                )
