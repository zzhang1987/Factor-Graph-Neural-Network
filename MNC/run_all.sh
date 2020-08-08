set -x

make clean
rm ansi/*.o
rm radford/*.o
make s2t
make t2y
make y2b
make zb2x 
./s2t -sfile _s/48 -k 48 -n 48 -Gfile codes/96.3.963/G -smn 1 -tfile _t/96
./t2y -tfile _t/96 -yfile _y/96 -gcx 1.32 -seed 322457 -n 96
./y2b -yfile _y/96 -n 96 -bfile _b/96 -gcx 1.32
./zb2x -bfile _b/96 -zfixed 0 -k 48 -n 48 -Afile codes/96.3.963/A2 -xfile _x/48gc -xso 1 -bndloops 100
diff _s/48 _x/48gc
