import subprocess

# Create zdipot.in file with the setup configuration
filename = 'test2.ss'
bincl = 50.
bvsin = 68.2
bvrad = -0.8
beta = -0.0015 
gamma = 0.0055
outputname = 'zdipot_search5'
f = zdipot_in(filename, outputname, incl = '$1', vsin = '$2', vrad = '$3', beta=beta, gamma=gamma, save_syntetic_spec='f')

incl = bincl + N.linspace(-10, 10, 5)
vsin = bvsin + N.linspace(-0.2, 0.2, 5)
vrad = bvrad + N.linspace(-0.2, 0.2, 5)
with open('search5.txt','a') as file:
    file.write('# beta = %6.4f and gamma = %6.4f \n' % (beta, gamma))
    file.write('incl vsin vrad  chi2       s         sp_ph   test  cool    hot     chi2I      chi2V'+'\n')
    for iincl in incl:
        for ivsin in vsin:
            for ivrad in vrad:
                complete = subprocess.run(['sh','./'+outputname, '%2d' % iincl, '%4.1f' % ivsin, '%3.1f' % ivrad], stdout=subprocess.PIPE)
                lastline = complete.stdout.decode('utf-8').splitlines()[-6].replace('<','').replace('(','').replace(')','').split()
                data = rzdipot(lastline)
                file.write('%2d ' % iincl); file.write('%4.1f ' % ivsin); file.write('%3.1f ' % ivrad)
                [file.write(element+' ') for element in data]; file.write('\n')
