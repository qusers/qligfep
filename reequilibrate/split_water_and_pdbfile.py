waterfilelines = ["25 sphere\n"]
proteinfilelines = []

with open("../md_0000_1000.pdb") as pdbfile:
    for line in pdbfile:
        if list(filter(None,line.split(" ")))[0].strip("\n") == "TER":
            proteinfilelines.append(line)
        else:
            RES = list(filter(None,line.split(" ")))[3].strip("\n")
            if RES == "HOH":
                waterfilelines.append(line)
            else:
                proteinfilelines.append(line)

with open("water.pdb", "w") as waterfile:
    for line in waterfilelines:
        waterfile.write(line)

# print(waterfilelines)
with open("protein.pdb", "w") as waterfile:
    for line in proteinfilelines:
        waterfile.write(line)