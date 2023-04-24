from feature import Feature
from genbankparser import GenbankParser


if __name__ == "__main__":
    input_file = "C:/Users/a.jalali/Downloads/Programming 2/Final Assignment-programming 1/modules/CFTR_mRNA.gb"
    feature = Feature(input_file)
    gp = GenbankParser(feature)
    format = input("please enter u for uppercased or s for separate: ")
    while True:
        if format == "u":
            print(gp.extract_features_uppercased())
            break
        elif format == "s":
            print(gp.extract_features_separated())
            break
        else:
            format = input("Failure! Please enter u  or s : ")
