import string
import re


class A10(object):
    def __init__(self, fn=None):
        self.created = None
        self.project = None
        self.stationname = None
        self.lat = None
        self.long = None
        self.elev = None
        self.setupht = None
        self.transferht = None
        self.actualht = None
        self.gradient = None
        self.nominalAP = None
        self.polarx = None
        self.polary = None
        self.dffile = None
        self.olfile = None
        self.clock = None
        self.blue = None
        self.red = None
        self.date = None
        self.time = None
        self.timeoffset = None
        self.gravity = None
        self.setscatter = None
        self.precision = None
        self.uncertainty = None
        self.collected = None
        self.processed = None
        self.transferhtcorr = None
        self.comments = None
        if fn:
            self.read_project_dot_txt(fn)

    def read_project_dot_txt(self, filename):
        dtf = False
        olf = False
        skip_grad = False
        in_comments = 0
        project_file = open(filename, 'r', encoding='unicode_escape')
        data_array = []  # ['a']*32
        # Look for these words in the g file
        tags = re.compile(r'Project|Name|Created|Setup' +
                          r'|Transfer|Actual|Date|Time|TimeOffset|Nominal|Red' +
                          r'|Blue|Scatter|SetsColl|SetsProc|Precision|Total_unc')
        # 'Lat' is special because there are three data on the same line:
        # (Lat, Long, Elev)
        lat_tag = re.compile(r'Lat')

        # 'Polar' is also special, for the same reason
        pol_tag = re.compile(r'Polar')

        version_tag = re.compile(r'Version')

        # Apparently using a delta file is optional, it's not always written to the .project file
        delta_tag = re.compile(r'DFFile')
        ol_tag = re.compile(r'OLFile')
        rub_tag = re.compile(r'RubFrequency')
        grav_tag = re.compile(r'Grv')
        grad_tag = re.compile(r'Gradient')

        # This one, because "Gradient:" is repeated exactly in this section
        unc_tag = re.compile(r'Uncertainties')

        # This deals with multi-line comments
        comment_tag = re.compile(r'Comments')

        for line in project_file:
            # Change up some text in the g file to make it easier to parse
            # (remove duplicates, etc.)
            line = line.strip()
            line = line.replace('\n\n', '\n')
            line = line.replace(":  ", ": ")
            # Repeat to take care of ":   " (three spaces)
            line = line.replace(":  ", ": ")
            line = line.replace(":  ", ": ")
            line = line.replace("g Acquisition Version", "Acq")
            line = line.replace("g Processing ", "")
            line = line.replace("Project Name:", "Project")
            line = line.replace("File Created:", "Created")
            line = line.replace('Gravity Corrections', 'grvcorr')
            line = line.replace(" Height:", ":")
            line = line.replace("Delta Factor Filename:", "DFFile")
            line = line.replace("Ocean Load ON, Filename:", "OLFile")
            line = line.replace("Nominal Air Pressure:", "Nominal")
            line = line.replace("Barometric Admittance Factor:", "Admittance")
            line = line.replace(" Motion Coord:", "")
            line = line.replace("Set Scatter:", "Scatter")
            line = line.replace("Offset:", "ofst")
            line = line.replace("Time Offset (D h:m:s):", "TimeOffset")
            line = line.replace("Ocean Load:", "OLC")
            line = line.replace("Rubidium Frequency:", "RubFrequency")
            line = line.replace("Blue Lock:", "Blue")
            line = line.replace("Red Lock:", "Red")
            line = line.replace("Red/Blue Separation:", "Separation")
            line = line.replace("Red/Blue Interval:", "Interval")
            line = line.replace("Gravity Corrections", "Corrections")
            line = line.replace("Gravity:", "Grv:")
            line = line.replace("Number of Sets Collected:", "SetsColl")
            line = line.replace("Number of Sets Processed:", "SetsProc")
            line = line.replace("Polar Motion:", "PolMotC")  # This is the PM error, not the values
            line = line.replace("Barometric Pressure:", "")
            line = line.replace("System Setup:", "")
            line = line.replace("Total Uncertainty:", "Total_unc")
            line = line.replace("Measurement Precision:", "Precision")
            line = line.replace(":", "", 1)
            line = line.replace(",", "")
            line_elements = line.split(" ")

            # Look for tags
            tags_found = re.search(tags, line)
            lat_tag_found = re.search(lat_tag, line)
            pol_tag_found = re.search(pol_tag, line)
            comment_tag_found = re.search(comment_tag, line)
            version_tag_found = re.search(version_tag, line)
            delta_tag_found = re.search(delta_tag, line)
            ol_tag_found = re.search(ol_tag, line)
            grav_tag_found = re.search(grav_tag, line)
            unc_tag_found = re.search(unc_tag, line)
            grad_tag_found = re.search(grad_tag, line)
            rub_tag_found = re.search(rub_tag, line)

            if unc_tag_found is not None:
                skip_grad = True

            if grad_tag_found is not None:
                if not skip_grad:
                    data_array.append(line_elements[1])

            # Old g versions don't output Time Offset, which comes right before gravity
            if grav_tag_found is not None:
                if version < 5:
                    data_array.append('-999')
                data_array.append(line_elements[1])

            if delta_tag_found is not None:
                dtf = True
                df = " ".join(line_elements[1:])

            if ol_tag_found is not None:
                olf = True
                of = " ".join(line_elements[1:])

            if rub_tag_found is not None:
                if dtf:
                    data_array.append(df)
                else:
                    data_array.append('-999')
                if olf:
                    data_array.append(of)
                else:
                    data_array.append('-999')
                data_array.append(line_elements[1])

            if version_tag_found is not None:
                version = float(line_elements[1])

            if tags_found is not None:
                try:
                    data_array.append(line_elements[1])
                except:
                    data_array.append('-999')

            if lat_tag_found is not None:
                data_array.append(line_elements[1])
                data_array.append(line_elements[3])
                data_array.append(line_elements[5])
                # This accomodates old versions of g. If these data are to be published,
                # though, they should be reprocessed in a more recent version.
                if version < 5:
                    data_array.append('-999')  # Setup Height
                    data_array.append('-999')  # Transfer Height
                    data_array.append('-999')  # Actual Height

            if pol_tag_found is not None:
                data_array.append(line_elements[1])
                data_array.append(line_elements[3])

            if in_comments > 0:
                comments += line
                if in_comments > 1:
                    comments += ' | '
                in_comments += in_comments

            if comment_tag_found is not None:
                in_comments = 1
                comments = ''

        data_array.append(comments)

        # Old g versions don't output transfer height correction
        if version < 5:
            data_array.append('-999')
        project_file.close()

        self.created = data_array[0]
        self.project = data_array[1]
        self.stationname = data_array[2]
        self.lat = data_array[3]
        self.long = data_array[4]
        self.elev = data_array[5]
        self.setupht = data_array[6]
        self.transferht = data_array[7]
        self.actualht = data_array[8]
        self.gradient = data_array[9]
        self.nominalAP = data_array[10]
        self.polarx = data_array[11]
        self.polary = data_array[12]
        self.dffile = data_array[13]
        self.olfile = data_array[14]
        self.clock = data_array[15]
        self.blue = data_array[16]
        self.red = data_array[17]
        date_elems = data_array[18].split('/')
        self.date = str(int(date_elems[2]) + 2000) + '-' + date_elems[0] + '-' + date_elems[1]
        self.time = data_array[19]
        self.timeoffset = data_array[20]
        self.gravity = data_array[21]
        self.setscatter = data_array[22]
        self.precision = data_array[23]
        self.uncertainty = data_array[24]
        self.collected = data_array[25]
        self.processed = data_array[26]
        self.transferhtcorr = data_array[27]
        self.comments = data_array[28]
