N=32*32
Bin(x)=N!/(x!*((N-x)!))

array A[10]

do for [i=1:10] {A[i]=i}

set xrange [0:10]
set sample 11
print Bin(2)
plot [1:10] '+' u ($0):(Bin($0)) t 'the points'
pause -1
