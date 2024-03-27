import os

directories = os.listdir('.')

# inputs = []

# for dirr in directories:
#     if dirr[:2] == 'md':
#         dirrr = dirr.strip('inp') + 'en'
#         inputs.append(dirrr)

# for i in reversed(inputs):
    # print(i)
    # print('time mpirun -np 16 $qdyn ' + i + ' > ' + i.replace('.inp','.log'))

for dirr in directories:
    # print(dirr)
    if dirr[0:7] == "mdfiles":
        for mdfiles in os.listdir(dirr):
            # print(mdfiles)

            if mdfiles[:2] == 'md' or mdfiles[:2] == "eq":
                newfilelines = []
                with open(dirr + "/" + mdfiles) as mdfile:
                    for line in mdfile:
                        # newfilelines.append(line)
                        if '[sequence_restraints]' in line:
                if 'shell_radius' in line:
                    newfilelines.append("shell_radius              25\n\n")
                else:
                    newfilelines.append(line)
                            break
                        else:
                            newfilelines.append(line)

                    with open(dirr + "/" + mdfiles, 'w') as new:
                        for line in newfilelines:
                            new.write(line)
