/*

	Smart memory allocator 

	Evolved over several years, starting with the initial SMARTALLOC
	code  for  AutoSketch in 1986, guided by the Blind Watchbreaker,
	John Walker.  Isolated in this general-purpose form in September
	of 1989.

*/

#include <stdio.h> 
#include <string.h>
#include <assert.h>

/*LINTLIBRARY*/
#ifdef lint
#include <memory.h>
extern char *sprintf();
#endif

#ifdef SMARTALLOC

typedef enum {False = 0, True = 1} Boolean;
typedef unsigned short ushort;

#define EOS      '\0'              /* End of string sentinel */
#define V	 (void)
#define min(a,b) ((a)<(b) ? (a) : (b))

/*  Queue data structures  */

/*  General purpose queue  */

struct queue {
	struct queue *qnext,	   /* Next item in queue */
		     *qprev;	   /* Previous item in queue */
};

/*  Memory allocation control structures and storage.  */

struct abufhead {
	struct queue abq;	   /* Links on allocated queue */
	unsigned ablen; 	   /* Buffer length in bytes */
	char *abfname;		   /* File name pointer */
	ushort ablineno;	   /* Line number of allocation */ 
};

static struct queue abqueue = {    /* Allocated buffer queue */
	&abqueue, &abqueue
};

static Boolean bufimode = False;   /* Buffers not tracked when True */

/*  Queue functions  */

void qinsert();
struct queue *qdchain();


/*  SMALLOC  --  Allocate buffer, enqueing on the orphaned buffer
		 tracking list.  */

static char *smalloc(fname, lineno, nbytes)
  char *fname;
  int lineno;
  unsigned nbytes;
{
	char *buf;
extern char *malloc();

	/* Note:  Unix	MALLOC	actually  permits  a zero length to be
	   passed and allocates a valid block with  zero  user	bytes.
	   Such  a  block  can	later  be expanded with realloc().  We
           disallow this based on the belief that it's better to  make
	   a  special case and allocate one byte in the rare case this
	   is desired than to miss all the erroneous occurrences where
	   buffer length calculation code results in a zero.  */

	assert(nbytes > 0);

	nbytes += sizeof(struct abufhead) + 1;
	if ((buf = malloc(nbytes)) != NULL) {
	   /* Enqueue buffer on allocated list */
	   qinsert(&abqueue, (struct queue *) buf);
	   ((struct abufhead *) buf)->ablen = nbytes;
	   ((struct abufhead *) buf)->abfname = bufimode ? NULL : fname;
	   ((struct abufhead *) buf)->ablineno = lineno;
	   /* Emplace end-clobber detector at end of buffer */
	   buf[nbytes - 1] = (((long) buf) & 0xFF) ^ 0xC5;
	   buf += sizeof(struct abufhead);  /* Increment to user data start */
	}
	return buf;
}

/*  SM_FREE  --  Update free pool availability.  FREE is never called
		 except  through  this interface or by actuallyfree().
		 free(x)  is  defined  to  generate  a	call  to  this
		 routine.  */

void sm_free(cp)
  char *cp;
{
	struct queue *qp;

	assert(cp != NULL);	   /* Better not release a null buffer, guy! */

	cp -= sizeof(struct abufhead);
	qp = (struct queue *) cp;

	/* The following assertions will catch virtually every release
           of an address which isn't an allocated buffer. */

	assert(qp->qnext->qprev == qp);   /* Validate queue links */
	assert(qp->qprev->qnext == qp);

	/* The following assertion detects storing off the  end  of  the
	   allocated  space in the buffer by comparing the end of buffer
	   checksum with the address of the buffer.  */

	assert(((unsigned char *) cp)[((struct abufhead *) cp)->ablen - 1] ==
	   ((((long) cp) & 0xFF) ^ 0xC5));

	V qdchain(qp);

	/* Now we wipe the contents of	the  just-released  buffer  with
           "designer  garbage"  (Duff  Kurland's  phrase) of alternating
	   bits.  This is intended to ruin the day for any miscreant who
           attempts to access data through a pointer into storage that's
	   been previously released. */

	V memset(cp, 0xAA, (int) ((struct abufhead *) cp)->ablen);

	free(cp);
}

/*  SM_MALLOC  --  Allocate buffer.  NULL is returned if no memory
		   was available.  */

char *sm_malloc(fname, lineno, nbytes)
  char *fname;
  int lineno;
  unsigned nbytes;
{
	char *buf;

	if ((buf = smalloc(fname, lineno, nbytes)) != NULL) {

	   /* To catch sloppy code that assumes  buffers  obtained  from
	      malloc()	are  zeroed,  we  preset  the buffer contents to
              "designer garbage" consisting of alternating bits.  */

	   V memset(buf, 0x55, (int) nbytes);
	}
	return buf;
}

/*  SM_CALLOC  --  Allocate an array and clear it to zero.  */

char *sm_calloc(fname, lineno, nelem, elsize)
  char *fname;
  int lineno;
  unsigned nelem, elsize;
{
	char *buf;

	if ((buf = smalloc(fname, lineno, nelem * elsize)) != NULL) {
	   V memset(buf, 0, (int) (nelem * elsize));
	}
	return buf;
}

/*  SM_REALLOC	--  Adjust the size of a  previously  allocated  buffer.
                    Note  that  the trick of "resurrecting" a previously
		    freed buffer with realloc() is NOT supported by this
		    function.	Further, because of the need to maintain
		    our control storage, SM_REALLOC must always allocate
		    a  new  block  and	copy  the data in the old block.
		    This may result in programs which make heavy use  of
		    realloc() running much slower than normally.  */

char *sm_realloc(fname, lineno, ptr, size)
  char *fname;
  int lineno;
  char *ptr;
  unsigned size;
{
	unsigned osize;
	char *buf;

	assert(size > 0);

	/*  If	the  old  block  pointer  is  NULL, treat realloc() as a
	   malloc().  SVID is silent  on  this,  but  many  C  libraries
	   permit this.  */

	if (ptr == NULL)
	   return sm_malloc(fname, lineno, size);

	/* If the old and new sizes are the same, be a nice guy and just
	   return the buffer passed in.  */

	osize = ((struct abufhead *) ptr)->ablen -
		(sizeof(struct abufhead) + 1);
	if (size == osize) {
	   return ptr;
	}

	/* Sizes differ.  Allocate a new buffer of the	requested  size.
           If  we  can't  obtain  such a buffer, act as defined in SVID:
	   return NULL from  realloc()	and  leave  the  buffer  in  PTR
	   intact.  */

	if ((buf = smalloc(fname, lineno, size)) != NULL) {
	   V memcpy(buf, ptr, (int) min(size, osize));
	   /* If the new buffer is larger than the old, fill the balance
              of it with "designer garbage". */
	   if (size > osize) {
	      V memset(buf + osize, 0x55, (int) (size - osize));
	   }

	   /* All done.  Free and dechain the original buffer. */

	   sm_free(ptr);
	}
	return buf;
}

/*  ACTUALLYMALLOC  --	Call the system malloc() function to obtain
			storage which will eventually be released
			by system or library routines not compiled
			using SMARTALLOC.  */

char *actuallymalloc(size)
  unsigned size;
{
extern char *malloc();

	return malloc(size);
}

/*  ACTUALLYCALLOC  --	Call the system calloc() function to obtain
			storage which will eventually be released
			by system or library routines not compiled
			using SMARTALLOC.  */

char *actuallycalloc(nelem, elsize)
  unsigned nelem, elsize;
{
extern char *calloc();

	return calloc(nelem, elsize);
}

/*  ACTUALLYREALLOC  --  Call the system realloc() function to obtain
			 storage which will eventually be released
			 by system or library routines not compiled
			 using SMARTALLOC.  */

char *actuallyrealloc(ptr, size)
  char *ptr;
  unsigned size;
{
extern char *realloc();

	return realloc(ptr, size);
}

/*  ACTUALLYFREE  --  Interface to system free() function to release
		      buffers allocated by low-level routines. */

void actuallyfree(cp)
  char *cp;
{
	free(cp);
}

/*  SM_DUMP  --  Print orphaned buffers (and dump them if BUFDUMP is
		 True). */

void sm_dump(bufdump)
  Boolean bufdump;
{
	struct abufhead *ap = (struct abufhead *) abqueue.qnext;

	while (ap != (struct abufhead *) &abqueue) {

	   if ((ap == NULL) ||
	       (ap->abq.qnext->qprev != (struct queue *) ap) || 
	       (ap->abq.qprev->qnext != (struct queue *) ap)) {
	      V fprintf(stderr,
                 "\nOrphaned buffers exist.  Dump terminated following\n");
	      V fprintf(stderr,
                 "  discovery of bad links in chain of orphaned buffers.\n");
	      V fprintf(stderr,
                 "  Buffer address with bad links: %lx\n", (long) ap);
	      break;
	   }

	   if (ap->abfname != NULL) {
	      unsigned memsize = ap->ablen - (sizeof(struct abufhead) + 1);
	      char errmsg[80];

	      V sprintf(errmsg,
                "Orphaned buffer:  %6u bytes allocated at line %d of %s\n",
		 memsize, ap->ablineno, ap->abfname
	      );
              V fprintf(stderr, "%s", errmsg);
	      if (bufdump) {
		 unsigned llen = 0;
		 char *cp = ((char *) ap) + sizeof(struct abufhead);

		 errmsg[0] = EOS;
		 while (memsize) {
		    if (llen >= 16) {
                       V strcat(errmsg, "\n");
		       llen = 0;
                       V fprintf(stderr, "%s", errmsg);
		       errmsg[0] = EOS;
		    }
                    V sprintf(errmsg + strlen(errmsg), " %02X",
		       (*cp++) & 0xFF);
		    llen++;
		    memsize--;
		 }
                 V fprintf(stderr, "%s\n", errmsg);
	      }
	   }
	   ap = (struct abufhead *) ap->abq.qnext;
	}
}

/*  SM_STATIC  --  Orphaned buffer detection can be disabled  (for  such
		   items  as buffers allocated during initialisation) by
		   calling   sm_static(1).    Normal   orphaned   buffer
		   detection  can be re-enabled with sm_static(0).  Note
		   that all the other safeguards still apply to  buffers
		   allocated  when  sm_static(1)  mode is in effect.  */

void sm_static(mode)
  int mode;
{
	bufimode = (Boolean) (mode != 0);
}

/*  Queue manipulation functions.  */


/*  QINSERT  --  Insert object at end of queue	*/

static void qinsert(qhead, object)
  struct queue *qhead, *object;
{
	assert(qhead->qprev->qnext == qhead);
	assert(qhead->qnext->qprev == qhead);

	object->qnext = qhead;
	object->qprev = qhead->qprev;
	qhead->qprev = object;
	object->qprev->qnext = object;
}


/*  QREMOVE  --  Remove object from queue.  Returns NULL if queue empty  */

static struct queue *qremove(qhead)
  struct queue *qhead;
{
	struct queue *object;

	assert(qhead->qprev->qnext == qhead);
	assert(qhead->qnext->qprev == qhead);

	if ((object = qhead->qnext) == qhead)
	   return NULL;
	qhead->qnext = object->qnext;
	object->qnext->qprev = qhead;
	return object;
}


/*  QDCHAIN  --  Dequeue an item from the middle of a queue.  Passed
		 the queue item, returns the (now dechained) queue item. */

static struct queue *qdchain(qitem)
  struct queue *qitem;
{
	assert(qitem->qprev->qnext == qitem);
	assert(qitem->qnext->qprev == qitem);

	return qremove(qitem->qprev);
}
#endif
