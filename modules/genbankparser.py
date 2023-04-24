from feature import Feature
import re
import math


class GenbankParser:
    """ This class gets an object of Feature class and parse it.
    It has two main methods for extracting freatures in uppercased or
    seperated format. According to the format, different function will be
    called. a method that validate the segment of a feature_title.
    """
    def __init__(self, feature: Feature):
        """ This method is called when GenbankParser object is created.
        the attributes are defined here. output file will be created too.
        """
        self.feature = feature      # <-- dependency is injected
        self.origin = self.feature.get_origin_sequence()
        self.input_file = self.feature.input_file
        self.definition = self.feature.get_definition()
        self.info_list = self.feature.get_list_of_features_info()
        self.output_file = self.create_output_file(self.definition)

    def create_output_file(self, seq):
        """ This method creates the output file and name it according to the
        input file name. it gets a string to write in the first line.
        """
        start = self.input_file.rfind("CFTR_")
        stop = self.input_file.rfind(".")
        # path + name (DNA.mRNA/protein)+  "_features.txt"
        output_file = self.input_file[:start] + self.input_file[start+5 : stop] + "_features.txt"
        with open(output_file, "w") as OF:
            OF.write(seq + "\n")
        return output_file

    def complement(self, dna_seq):
        """ This method gets a sequence and converts it to
        the reverse and complement format
        """
        comp_bases = {"a":"t", "c":"g", "t":"a", "g":"c", "A": "T", "T": "A", "C": "G", "G": "C"}
        reverse_dna = dna_seq[::-1]
        x = []
        for i in reverse_dna:
            x.append(comp_bases[i])
        Complementary_seq = ''.join(x)
        return Complementary_seq

    def trim_sequence(self, seq):
        """ This method trims the sequence.
        it gets a string and if the length is more thans 60 characters,
        writes the first 60charcter in the output_file  and goes to the
        next line and does it until the length is less than 60
        """
        a = 60
        length = math.ceil(len(seq) / a)
        with open(self.output_file, "a") as EF:
            for i in range(length):
                start = i*a
                stop = a *(i+1)
                EF.write(seq[start: stop] + "\n")
            if len(seq) % a == 0:
                EF.write("\n")

    def validate_segment_of_sequence(self, length, *num):
        """ This method that checks whether the segment of a
        feature_title exists. If the start or stop of segment is
        greater than of the length of sequence, it raise a error.
        """
        for n in num:
            if n > length:
#                 raise Errors
                raise ValueError(f"{n} is greater than of the length of sequence")

    def extract_join_uppercased(self, location):
        """this method make a final sequence from main and rest;only main
        segments are uppercased, sequence between them are lowercased.
        """
        list_num = [int(num) for num in re.findall(r"\d+", location)]
        start_origin_seq = 0
        Final_Sequence = ""
        for i in range(0, len(list_num), 2):
            start = list_num[i]
            stop = list_num[i + 1]
            stop_origin_seq = start - 1
            # define main and rest seperately
            main_Sequence = self.origin[start-1:stop]
            rest_Sequence = self.origin [start_origin_seq :stop_origin_seq]
            main_Sequence = main_Sequence.upper()
            Final_Sequence = Final_Sequence + rest_Sequence + main_Sequence
            start_origin_seq = stop
        return Final_Sequence

    def extract_join_separated(self, location):
        """this method make a final sequence by joining all segments.
        """
        Final_Sequence = ''
        list_num = [int(num) for num in re.findall(r"\d+", location)]
        if len(list_num) % 2 == 0:
            for i in range(0, len(list_num), 2):
                start = list_num[i]
                stop = list_num[i + 1]
                Final_Sequence = Final_Sequence + self.origin[start-1:stop]
        return Final_Sequence

    def extract_normal_uppercased(self, location):
        """this method make a final sequence by making the segment uppercased
        and adding the rest of previous sequence to it.
        """
        try:
            list_num = [int(num) for num in re.findall(r"\d+", location)]
            start = list_num[0]
            stop_origin_seq = start - 1
            if len(list_num) == 2:
                stop = list_num[1]
                # at first the location of segment is validated
                self.validate_segment_of_sequence(len(self.origin ), start, stop, stop_origin_seq)
                main_Sequence = self.origin[start-1:stop]
            else:
                self.validate_segment_of_sequence(len(self.origin), start)
                main_Sequence = self.origin[start-1]
            main_Sequence = main_Sequence.upper()
            rest_Sequence = self.origin[0:stop_origin_seq]
            Final_Sequence = rest_Sequence + main_Sequence
        except ValueError:
            Final_Sequence = ""
        return Final_Sequence

    def extract_normal_separated(self, location):
        """extract the segment from origin sequence.
        """
        try:
            list_num = [int(num) for num in re.findall(r"\d+", location)]
            start = list_num[0]
            if len(list_num) == 2:
                stop = list_num[1]
                # at first the location of segment is validated
                self.validate_segment_of_sequence(len(self.origin), start, stop)
                Final_Sequence = self.origin[start-1:stop]
            else:
                self.validate_segment_of_sequence(len(self.origin), start - 1)
                Final_Sequence = self.origin[start-1]
        # if start or stop are out of origin length, nothing is written in file
        except ValueError:
            Final_Sequence = ""
        return Final_Sequence

    def write_in_file(self, current_object, Final_Sequence):
        """ this method writes title of feature, related sequence of feature
        and its expaination(next line of numbers) in output file.
        """
        with open(self.output_file, "a") as EF:
            EF.write("\n>" + self.info_list[current_object] + ' ' + self.info_list[current_object+1][1] + "\n")
        self.trim_sequence(Final_Sequence)

    def extract_features_uppercased(self):
        """ This method will be called if the requested format is uppercased.
        It analyzes the items of info_list from self.get_list_of_features_info.
        extracts the numbers(locations). at first, finds the rest_Sequence and
        main_Sequence according to the location.
        """
        number_of_objects  = len(self.info_list)
        main_Sequence = ''
        current_object = 0
        while current_object < number_of_objects :
            main_Sequence = ''
            location = self.info_list[current_object+1][0]
            # if finds 'join' in the line, calls extract_join_uppercased;
            if 'join' in location:
                Final_Sequence = self.extract_join_uppercased(location)
            # if finds 'order', nothing should be extracted.
            elif 'order' in location:
                   Final_Sequence = ""
            # normal case
            else:
                Final_Sequence = self.extract_normal_uppercased(location)
            # if finds 'complement', calls complement method.
            if 'complement' in location:
                    Final_Sequence = self.complement(Final_Sequence)
            # write in ouput in a trimmed way
            if Final_Sequence != "":
                self.write_in_file(current_object, Final_Sequence)
            # go for next feature
            current_object += 2
        # if all features are not parsed returns an error
        try:
            if current_object == number_of_objects:
                return "GenBank Features extracted Successfully"
            else:
                raise TypeError("not correct file")
        except TypeError:
            return "GenBank Features is NOT extracted completely"

    def extract_features_separated(self):
        """ This method will be called if the requested format in separated.
        it analyzes the items of info_list from self.get_list_of_features_info.
        It analyzes the second items of self.list_of_feature_objects_info.
        extracts the numbers(locations). finds the Final_Sequence according
        to the location.
        """
        number_of_objects  = len(self.info_list)
        current_object = 0
        while current_object < number_of_objects :
            location = self.info_list[current_object+1][0]
            # if finds 'join' in the line, calls extract_join_separated;
            if 'join' in location:
                Final_Sequence = self.extract_join_separated(location)
            # if finds 'order', nothing should be extracted.
            elif 'order' in location:
                    Final_Sequence = ""
            # normal case
            else:
                Final_Sequence = self.extract_normal_separated(location)
            # if finds 'complement', calls complement method.
            if 'complement' in location:
                Final_Sequence = self.complement(Final_Sequence)
            # write in ouput in a trimmed way
            if Final_Sequence != "":
                self.write_in_file(current_object, Final_Sequence)
            # go for next feature
            current_object += 2
        # if all features are not parsed returns an error
        try:
            if current_object == number_of_objects:
                return "GenBank Features extracted Successfully"
            else:
                raise TypeError("not correct file")
        except TypeError:
            return "GenBank Features is NOT extracted completely"
