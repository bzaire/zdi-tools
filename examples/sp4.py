import subprocess

# Create zdipot.in file with the setup configuration
filename = 'test2.ss'
incl = 50.
vsin = 68.2
vrad = -0.8
outputname = 'zdipot_search4'
f = zdipot_in(filename, outputname, incl, vsin, vrad=vrad, beta='$1', gamma = '$2', save_syntetic_spec='f')

beta = N.linspace(-.0023, -.0003, 10)
gamma = N.linspace(0.0046, 0.0066, 10)
with open('search4.txt','a') as file:
#    file.write('beta gamma chi2       s         sp_ph   test  cool    hot     chi2I      chi2V'+'\n')
    for ib in beta:
        for ig in gamma:
            complete = subprocess.run(['sh','./'+outputname, '%6.4f' % ib, '%6.4f' % ig], stdout=subprocess.PIPE)
            lastline = complete.stdout.decode('utf-8').splitlines()[-6].replace('<','').replace('(','').replace(')','').split()
            data = rzdipot(lastline)
            file.write('%6.4f ' % ib); file.write('%6.4f ' % ig)
            [file.write(element+' ') for element in data]; file.write('\n')
