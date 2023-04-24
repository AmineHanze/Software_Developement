import re


class Feature:
    """ This class reads input file and extracts the necessary parts
    using different functions.
    """
    def __init__(self, input_file):
        """ initialize the attributes;This method is called
        when Feature object is created. input_file in the only
        attribute of Feature class when get instantiated.
        """
        self.input_file = input_file
        self.extract_info_from_input_file()

    def make_list_of_features(self, iterator):
        """creates a list of Feature objects and a list of their info:
        [Feature1, [sequence location , explanation], Feature2,[...], ...]
        """
        self.list_of_features_info = []
        line = next(iterator)
        while line.find('ORIGIN', 0, 20) == -1:
            line2 = line[:21].strip()
            if line2:
                self.list_of_features_info.append(line2)
                info = []
                info_seq = ''
                # when it finds "(", looks for ')' to save all locations.
                if '(' in line:
                    while line.find(')', 21, 80) == -1:
                        info_seq = info_seq + line[21:]
                        line = next(iterator)
                info_seq = info_seq + line[21:]
                info_explain = next(iterator).strip()
                # make this list:[sequence location, explanation]
                info.append(info_seq)
                info.append(info_explain)
                self.list_of_features_info.append(info)
            line = next(iterator)
        return iterator

    def make_origin_sequence(self, line, iterator):
        """ Extract a string of all lines in ORIGIN part
        """
        self.origin_sequence = ""
        while line != "//":
            # remove line numbers
            line = re.sub(r'[0-9]', '', line)
            # remove spaces
            line = re.sub('\s', '', line)
            self.origin_sequence = self.origin_sequence + line
            line = next(iterator).strip()
        # remove "iflqpe" letters from protein
        if "m" in self.origin_sequence:
            self.origin_sequence = re.sub(r'[iflqpe]', '', self.origin_sequence)

    def extract_info_from_input_file(self):
        """this method extracts all needed information from
        input file by openning the file in read mode just once.
        read it line by line and only save data in variables.
        """
        self.sequence = ""
        with open(self.input_file, 'r') as f:
            iterator = iter(f)
            line = next(iterator)
            # extracts data before last line of file that contains "//"
            while line.find('//', 0, 2) == -1:
                if line.find('DEFINITION', 0, 20) == 0:
                    self.definition = line.strip().replace("DEFINITION  ", "")
                    if not self.definition.endswith("."):
                        self.definition = self.definition + " " + next(iterator).strip()
                elif line.find('FEATURES', 0, 20) == 0:
                    iterator = self.make_list_of_features(iterator)
                elif line.find('        1', 0, 20) == 0:
                    self.make_origin_sequence(line, iterator)
                    break
                line = next(iterator)

    def get_origin_sequence(self):
        """this method returns the string of origin sequence
        """
        return self.origin_sequence

    def get_definition(self):
        """this method returns the string of definition
        """
        return self.definition

    def get_list_of_features_info(self):
        """this method returns the list of features info
        """
        return self.list_of_features_info

    def __iter__(self, file):
        """refere to __next__ function
        """
        return self.__next__()

    def __next__(self):
        """ return one line of file
        """
        try:
            yield f.readline()
        except StopIteration:
            return "no more lines"
