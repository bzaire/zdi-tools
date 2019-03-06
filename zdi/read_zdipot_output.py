def rzdipot(lastline):
    """ This funtion reads data format of zdipot routine.
        Parameters:
            lastline: str
        Return:
            list with the formated outputs of zdipot
    """
    chi2 = lastline[2]     # Chi square
    s = lastline[3]        # Entropy
    sp_ph = lastline[4]    # Spot coverage
    test = lastline[5]     # Test
    cool = lastline[7]     # Cool spot
    hot = lastline[9]     # Hot spot
    chi2I = lastline[14]   # individual chi square
    chi2V = lastline[15]   #
    return chi2, s, sp_ph, test, cool, hot, chi2I, chi2V
