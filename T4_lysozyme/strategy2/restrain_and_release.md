ligand in protein is initially restrained and then released near end point. Equilibration is performed at lamda 0.5 with charges off. Charges are added in a seperate FEP leg

We have tried 3 different alternatives in this scheme.

1: softcore is on in FEP1 and remains on

2: softcore is turned off during FEP1 (this does not work in Q I found out later)

3: no softcore is used at all. end point values are not taken into account in the calculation of deltaG

The uploaded inputfiles are the ones without use of softcore. Hoever they can easily be transformed to the other schemes.

![ strategy2_start_in_middle_cycle](strategy2_start_in_middle_cycle.jpg)