#!/usr/local/bin/perl
#
# code to read from special file and spit out pieces of code for 
# setting defaults and for printing the usage command and for 
# reading the command line
#
# the argument should be a file in which items have the form
# variable    type default flag  short_name extra_action long_name
# nc->train_n d    100     -nnn  n          -            number to train on
# nc->infile  s    -       -nnin nninfile   nc->read=2;  file to get weights from
#
# each line is turned into a structure entry   	(str)
#     	a default-setting entry			(def)
#     	a usage-printing entry			(usg)
# 	and a command line reading entry	(clr)
# of these, the S, U and C can be disabled by >No S, etc.
# and the D can be disabled by -.
#

$verbose = 1 ; 
$comcol  = 30 ; 
$longcol = 30;
$structcol = 20;
$comlength = 18 ; 
$restlength = 20 ; 
$adfl = 6 ; 

eval "\$$1=\$2" while $ARGV[0]=~ /^(\w+)=(.*)/ && shift;

$file = $ARGV[0] ;
$filedef = $ARGV[0]."_def.c" ;
$fileusg = $ARGV[0]."_usg.c" ;
$fileclr = $ARGV[0]."_clr.c" ;
$filestr = $ARGV[0]."_str.h" ;
$charnum = "[100]" ;

%type_to_print = ( 'd' , 'd' ,
		   'f' , '9.3g' ,
		   'ld' , 'ld' ,
		  's' , 's' ) ;
%type_to_scan = ( 'd' , 'd' ,
		 'f' , 'lf' ,
		 'ld' , 'ld' ,
		 's' , 's' ) ;
%type_to_struct = ( 'd' , 'int' ,
		   'f' , 'double' ,
		   'ld' , 'long int' ,
		   's' , 'char' ) ;

open ( DEF , "> $filedef" ) ;
open ( USG , "> $fileusg" ) ;
open ( CLR , "> $fileclr" ) ;
open ( STR , "> $filestr" ) ;

$doC = 1 ; $doU = 1 ; $doS = 1 ; 

while ( <> ) {			# 
    if ( /^\>/ ) {		# 
				# Command to interpret 
	if ( /^\>P/ ) { # print something somewhere ; 
	    s/^\>P// ; 
	    if ( s/^C// ) { print CLR ; }
	    if ( s/^U// ) { print USG ; }
	    if ( s/^S// ) { print STR ; } 
	    if ( s/^D// ) { print DEF ; } 
	} else {
	    if ( /[Nn][Oo]/ ) { $set = 0 ; } else { $set = 1 ; }
	    if ( /C/ ) { $doC = $set ; }
	    if ( /U/ ) { $doU = $set ; }
	    if ( /S/ ) { $doS = $set ; } # 
	    if (verbose>=2 ) { print STDERR $doC,$doU,$doS,"\n"; }
	}
    } elsif ( !/^\#|^%|^\/\*/ ) {

	($var , $type , $def , $flag , $short , $action , @resta) = split ;

	$rest = join ( " " , @resta ) ; 
	$rest = $rest." " x ( $restlength - length($rest) ) ; 
	if ( $def ne "-" ) {	# 
	    $tmps = "  $var = $def ;" ;
	    $fill = " " x ( $comcol - length($tmps) ) ; 
	    $tmpc = "$flag $short" ; 
	    $fillc = " " x ( $comlength - length($tmpc) ) ; 
	    print DEF "$tmps$fill/* $tmpc$fillc */\n" ;
	    $usg = $type_to_print{$type} ; 
	    $body = "$flag $short <\%$usg>" ; 
	    $tail = ", $var" ;
	    $extra = ($type eq "f") ? 7 : ( length($def) - 2 ) ; 
	} else {
	    $body = "$flag $short" ; # 
	    $tail = "" ;
	    $extra = 0 ; 
	}			# 
	$fill = " " x ( $longcol - length($body) - $extra ) ; 
	$qbody = '"'.$body.$fill.'('.$rest.')"' ;

	if ( $doU ) { print USG "  DNT; fprintf( fp, $qbody$tail);\n" ; }
	
	$qflag = '"'.$flag.'"' ; 
	
	if ( $doC ) {
	    print CLR "    else if ( strcmp (argv[i], $qflag) == 0 ) {\n" ;
	    if ( $type eq "s" ) {
		print CLR "       if ( i + 1 == argc ) { ERROR1; }\n" ; 
		print CLR "       else {\n" ; 
		print CLR "         strcpy($var, argv[++i]);\n" ;
		if ( $action ne "-" ) {
		    print CLR "         $action\n" ;
		}			# 
		print CLR "       }\n" ; 
	    } else {
		$qclr =  '"%'.$type_to_scan{$type}.'"' ; 
		print CLR "      if ( i + 1 == argc ) { ERROR1; }\n" ;
		print CLR "      else {\n" ; 
		print CLR "        cs *= sscanf(argv[++i], $qclr, &($var));\n" ;
		if ( $action ne "-" ) {
		    print CLR "         $action\n" ;
		}			# 
		print CLR "       }\n" ; 
	    }
	    print CLR "    }\n" ; 
	}
	
	$devar = $var ; 
	$devar =~ s/^.*\-\>// ;
	if ( $type eq "s" ) {
	    $body = "  ".$type_to_struct{$type}." ".$devar.$charnum." ;";
	} else {
	    $body = "  ".$type_to_struct{$type}." ".$devar." ;";
	}
	$fill = " " x ( $structcol - length($body) ) ; 
	$angledef = "<$def>" ;
	$angledef .= " " x ( $adfl - length ( $angledef ) ) ; 
	if ( $doS ) { print STR "$body$fill/* $rest $angledef */\n" ;}
    }
}

