#include "../ansi/r.h"
#include "../ansi/rand.h"
#include "../ansi/mynr.h"

void main(int argc, char *argv[])
{
  int i;
  long int seed = 123 ;

  ran_seed(seed) ;

  for ( i = 1 ; i <= 10000 ; i++ ) {
    printf("%g %g\n",rann(), ranc() );
  }
}
 
/*

histo.p min=-4 max=4 c=1 _out > _his
histo.p min=-4 max=4 c=2 _out > _hist
gnuplot
plot "_his" u 1:2 w linesp, "_his" u 1:3 w l
replot '_hist' u 1:2 w linesp, '_hist' u 1:3 w l

*/
