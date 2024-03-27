import os

directories = os.listdir('.')

# inputs = []

# for dirr in directories:
#     if dirr[:2] == 'md':
#         # dirrr = dirr.strip('inp') + 'en'
#         inputs.append(dirr)

# for i in inputs:
#     # print(i)
#     print('time mpirun -np 24 $qdyn ' + i + ' > ' + i.replace('.inp','.log'))

for dirr in directories:
    start = 0
    # print(dirr)
    if dirr[:2] == 'md' or dirr[:2] == "eq":
        newfilelines = []
        with open(dirr) as mdfile:
            for line in mdfile:
                # newfilelines.append(line)
                if '[sequence_restraints]' in line:
                    if dirr[:3] == "eq1" or dirr[:3] == "eq2":
                        newfilelines.append(line)
                        newfilelines.append('2597   2620   10.0 0  0\n\n')
                        newfilelines.append('[distance_restraints]\n')
                    else:
                        newfilelines.append(line)
                        newfilelines.append('2597   2620    5.0 0  0\n\n')
                        newfilelines.append('[distance_restraints]\n')
                    break
                else:
                    newfilelines.append(line)

            with open(dirr, 'w') as new:
                for line in newfilelines:
                    new.write(line)