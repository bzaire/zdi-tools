import subprocess

# Create zdipot.in file with the setup configuration
filename = 'test2.ss'
incl = 80.
vsin = 80.2
f = zdipot_in(filename, incl, vsin, vrad='$1')

vrad = N.linspace(-0.6, 1, 15)
#vrad = -0.8
with open('search2.txt','a') as file:
    file.write('incl vsin vrad chi2       s         sp_ph   test  cool    hot     chi2I      chi2V'+'\n')
    for ivrad in vrad:
        complete = subprocess.run(['sh','./zdipot_'+filename[:-3], str(ivrad)], stdout=subprocess.PIPE)
        lastline = complete.stdout.decode('utf-8').splitlines()[-6].replace('<','').replace('(','').replace(')','').split()
        data = rzdipot(lastline)
        file.write(str(incl)+' '); file.write(str(vsin)+' '); file.write('%03.1f ' % ivrad)
        [file.write(element+' ') for element in data]; file.write('\n')
