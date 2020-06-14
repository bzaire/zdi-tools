import subprocess

# Create zdipot.in file with the setup configuration
filename = 'test2.ss'
incl = 80.
vsin = 80.2
outputname = 'zdipot_search3'
f = zdipot_in(filename, outputname, incl, vsin, beta='$1', gamma = '$2')

beta = N.linspace(-.0080, -.0005, 15)
gamma = N.linspace(0.0005, 0.0070, 20)
with open('search3.txt','a') as file:
    file.write('# incl = %2d and vsin = %4.1f' % (incl, vsin))
    file.write('# beta gamma chi2       s         sp_ph   test  cool    hot     chi2I      chi2V'+'\n')
    for ib in beta:
        for ig in gamma:
            complete = subprocess.run(['sh','./'+outputname, '%6.4f' % ib, '%6.4f' % ig], stdout=subprocess.PIPE)
            lastline = complete.stdout.decode('utf-8').splitlines()[-6].replace('<','').replace('(','').replace(')','').split()
            data = rzdipot(lastline)
            file.write('%6.4f ' % ib); file.write('%6.4f ' % ig)
            [file.write(element+' ') for element in data]; file.write('\n')
