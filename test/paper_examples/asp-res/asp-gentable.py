from numpy import *

def frexp10(x):
    exp = int(log10(x))
    return x/10**exp,exp

def w_f_str(x):
    s = "{:.2f}".format(x)
    return s 

def w_g_str(x):
    s = '{:g}'.format(float('{:.3g}'.format(x)))
    return s

def w_g_str2(x):
    s = '{:g}'.format(float('{:.4g}'.format(x)))
    return s


norder = 4
z = loadtxt('stell_asp_mac_order'+str(norder)+'_iprec2.txt')
inc_header = 0
inc_end = 0

[m,n] = shape(z)

row_headers = ['\\aavg','\\Npat','m','\\alpha','\\snear','\\slp']

icols = [1,2,5,7,9,11]
data = z[:,icols]
dirname = ''


f = open(dirname+'asp-table-res'+str(norder)+'.tex','w')
if(inc_header == 1):
    f.writelines('\\documentclass{article}\n')
    f.writelines('\\usepackage{tikz}\n')
    f.writelines('\\usepackage[skip=-0.5\\baselineskip]{subcaption}\n')

    f.writelines('\\begin{document}\n')
    
f.writelines('\\begin{table}\n')
f.writelines('\\begin{center}\n')
f.writelines('\\[\\begin{array}{|c|c|c|c|c|c|}\n')
f.writelines('\\hline\n')
for i in range(size(icols)):
    s1 = row_headers[i]
    if(i != 1):
        for j in range(m):
            s1 = s1+' & '+ w_g_str(data[j,i])
        s1 = s1+'\\\ \hline \n'
    if(i == 1):
        for j in range(m):
            s1 = s1+' & '+ w_g_str2(data[j,i])
        s1 = s1+'\\\ \hline \n'
    f.writelines(s1)

f.writelines('\\end{array}\\]\n')
f.writelines('\\end{center}\n')
f.writelines('\\caption{Performance as a function of average aspect ratio $\\aavg$}\n')
f.writelines('\\label{tab:numerical-asp-res}\n')
f.writelines('\\end{table}\n')
if(inc_end == 1):
    f.writelines('\\end{document}\n')
f.close()
#
