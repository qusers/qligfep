import os

directories = os.listdir('.')

# inputs = []

# for dirr in directories:
#     if dirr[:2] == 'md':
#         dirrr = dirr.strip('inp') + 'en'
#         inputs.append(dirrr)

foldername = "large_systematic_changes/T4_1_1_1_0_1"

# for i in reversed(inputs):
    # print(i)
    # print('time mpirun -np 16 $qdyn ' + i + ' > ' + i.replace('.inp','.log'))
for dirrr in directories:
    if os.path.isdir(dirrr):
        directories2 = os.listdir(dirrr + "/inputfiles")
        for file in directories2:
            # print(file)
            if file == 'runSNELLIUS_only_analysis.sh':
                # print('yes')
                newfilelines = []
                with open(dirrr + "/inputfiles/" + file) as mdfile:
                    for mdfile_lines in mdfile:
                        if "python" in mdfile_lines:
                            # print(mdfile_lines)
                            # print(mdfile_lines[:-2] + "1\n")
                            newfilelines.append(mdfile_lines.replace("python", "~/anaconda3/envs/QligFEP/bin/python"))
                        else:
                            newfilelines.append(mdfile_lines)

                with open(dirrr + "/inputfiles/runSNELLIUS_only_analysis.sh", 'w') as new:
                    for line in newfilelines:
                        new.write(line)