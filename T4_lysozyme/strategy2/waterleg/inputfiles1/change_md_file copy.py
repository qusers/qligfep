import os

directories = os.listdir('.')

# inputs = []

# for dirr in directories:
#     if dirr[:2] == 'md':
#         dirrr = dirr.strip('inp') + 'en'
#         inputs.append(dirrr)

# for i in reversed(inputs):
#     print(i)
    # print('time mpirun -np 24 $qdyn ' + i + ' > ' + i.replace('.inp','.log'))

for dirr in directories:
    start = 0
    # print(dirr)
    if dirr[:2] == 'md':
        newfilelines = []
        with open(dirr) as mdfile:
            for line in mdfile:
                # newfilelines.append(line)
                if '[sequence_restraints]' in line:
                    newfilelines.append(line)
                    newfilelines.append('1      12       5.0 0 0\n\n')
                    newfilelines.append('[distance_restraints]\n')
                    break
                else:
                    newfilelines.append(line)

            with open(dirr, 'w') as new:
                for line in newfilelines:
                    new.write(line)