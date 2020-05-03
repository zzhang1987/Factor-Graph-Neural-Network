From sanger.ac.uk!rfs Wed Jul 13 15:59:08 1994
From: rfs@sanger.ac.uk
To: mackay@mrao.cam.ac.uk
Subject: immediate input

Here is the c program you requested - further info from man.termio and
man.ioctl

Obviously all the printf lines are not necessary.

/*  Last edited: Jul 13 15:46 1994 (rfs) */
#include <stdio.h>
#include <termio.h>
main()
{
  int c;

  struct termio termdat,termdatnew;
  struct termio *tdp;
  tdp=&termdat;

  c=ioctl((*stdin)._file,TCGETA,tdp);

  termdatnew=termdat;
  
  printf("ioctl output (should be zero) is %d\n\r",c); 

  printf("termdat.c_cflag is 0%o\n\r",termdat.c_cflag);
  printf("termdat.c_iflag is 0%o\n\r",termdat.c_iflag);
  printf("termdat.c_oflag is 0%o\n\r",termdat.c_oflag);
  printf("termdat.c_lflag is 0%o\n\r",termdat.c_lflag);
  printf("termdat.c_line is &%x\n\r",termdat.c_line);

  printf("termdat.c_cc[INTR] is &%x\n\r",termdat.c_cc[0]);
  printf("termdat.c_cc[QUIT] is &%x\n\r",termdat.c_cc[1]);
  printf("termdat.c_cc[ERASE] is &%x\n\r",termdat.c_cc[2]);
  printf("termdat.c_cc[KILL] is &%x\n\r",termdat.c_cc[3]);
  printf("termdat.c_cc[EOF] is &%x\n\r",termdat.c_cc[4]);
  printf("termdat.c_cc[EOL] is &%x\n\r",termdat.c_cc[5]);
  printf("termdat.c_cc[EOL2] is &%x\n\r",termdat.c_cc[6]);
  printf("termdat.c_cc[SWTCH] is &%x\n\r",termdat.c_cc[7]);
  printf("termdat.c_cc[START] is &%x\n\r",termdat.c_cc[8]);
  printf("termdat.c_cc[STOP] is &%x\n\r",termdat.c_cc[9]);
  printf("termdat.c_cc[SUSP] is &%x\n\r",termdat.c_cc[10]);
  printf("termdat.c_cc[REPRINT] is &%x\n\r",termdat.c_cc[12]);
  printf("termdat.c_cc[DISCARD] is &%x\n\r",termdat.c_cc[13]);
  printf("termdat.c_cc[WERASE] is &%x\n\r",termdat.c_cc[14]);
  printf("termdat.c_cc[LNEXT] is &%x\n\r",termdat.c_cc[15]);

  termdatnew.c_cc[VMIN]=1;
  termdatnew.c_cc[VTIME]=1;
  termdatnew.c_lflag&=~ICANON;
  termdatnew.c_lflag&=~ECHO;
  tdp=&termdatnew;
  c=ioctl((*stdin)._file,TCSETA,tdp); /* QQQ */
    
  printf("Hi, this program should now start in immediate mode \r\n");
  printf(" - to see what would have happened without fiddling, comment out the line containing QQQ \r\n");
  c=getchar();
  printf("char got was %x\n\x0D",c);    

  while ((c=getchar())!=4)
    putchar(c);

  printf("\n\rYour terminal has been LEFT in immediate mode!!!-but your shell may reset things\n\r");

/*add these line to reverse changes 
  tdp=&termdat;
  c=ioctl((*stdin)._file,TCSETA,tdp); 
*/

}










