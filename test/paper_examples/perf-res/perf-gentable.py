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

z = loadtxt('stell20_perf_mac.txt')


varname1 = '$p$'
varname2 = '$\\varepsilon$'
inc_header = 0
inc_end = 0

[m,n] = shape(z)

epsvals = zeros(4)
epsvals = z[0:4,1]

epstxt = ['5 \\cdot 10^{-2}','5 \\cdot 10^{-3}','5 \\cdot 10^{-6}','5 \\cdot 10^{-9}']

pvals = zeros(5)
pvals = z[0:m:4,0]

icols = [5,7,9,11]
dirname = ''
#icap = ['$m$','$\\alpha$','$\\snear$','$\\slp$']
#icap = ['$m$','$\\alpha$','$s_1$','$s_2$']
icap = ['$m$','$\\alpha$','$\\snear$','$\\slp$']


f = open(dirname+'perf-table-res.tex','w')
if(inc_header == 1):
    f.writelines('\\documentclass{article}\n')
    f.writelines('\\usepackage{tikz}\n')
    f.writelines('\\usepackage[skip=-0.5\\baselineskip]{subcaption}\n')

    f.writelines('\\begin{document}\n')
    
f.writelines('\\begin{table}\n')

for i in range(size(icols)):
    data = z[:,icols[i]]
    data = data.reshape(size(pvals),size(epsvals))

    f.writelines('\\begin{minipage}{0.5\\linewidth}\n')
    f.writelines('\\begin{center}\n') 
    f.writelines('\\[\\begin{array}{|c|c|c|c|c|}\n')
    f.writelines('\\hline')
    f.writelines('\\tikz{\n')
    f.writelines('\\node[below left, inner sep=1pt] (v1) {'+varname1+'};\n')
    f.writelines('\\node[above right, inner sep=1pt] (v2) {'+varname2+'};\n')
    f.writelines('\\draw(v1.north west|-v2.north west) -- (v1.south east-|v2.south east);}\n')

    s1 = ''
    for j in range(size(epsvals)):
        s1 = s1+' & '+ epstxt[j]
    s1 = s1+'\\\ \hline \n'
    f.writelines(s1)

    for j in range(size(pvals)):
        s1 = str(int(pvals[j])+1)
        for k in range(size(epsvals)):
            s1 = s1 + ' & ' + w_g_str(data[j,k])
        s1 = s1+ '\\\ \hline \n'
        f.writelines(s1)
    f.writelines('\\end{array}\\]\n')
    f.writelines('\\subcaption{'+icap[i]+'}\n')
    f.writelines('\\end{center}\n')
    f.writelines('\\end{minipage}\n')
    if(i==1):
        f.writelines('\\newline\n')
#f.writelines('\\caption{$m$,$\\snear$, and $\\slp$ as a function of $p$ and $\\varepsilon$}\n')
f.writelines('\\caption{$m$,$\\alpha$,$s_{1}$, and $s_{2}$ as a function of $p$ and $\\eta$}\n')
f.writelines('\\label{tab:numerical-perf-res}\n')
f.writelines('\\end{table}\n')
if(inc_end == 1):
    f.writelines('\\end{document}\n')
f.close()
#
