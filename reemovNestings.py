def reemovNestings(MiR):
    output=[];
    for i in MiR:
        if type(i) == list:
            reemovNestings(i)
        else:
            output.append(i)

    return output