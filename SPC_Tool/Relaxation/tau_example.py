
from mag_anal.Anal_Tau import CAL_RELAX


if __name__ == "__main__":

    """
    This script uses yaml as parameter input method
    """
    Input_filename = "inpfile/CFPs_v1.yml"
    """
    reading the input_file
    """

    TAU = CAL_RELAX(Input_filename)
    TAU.kernel()
    TAU.TP[0].Mod_Analyse_OB()