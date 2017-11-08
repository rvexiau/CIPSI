#include <time.h>  
#include <string.h>

 timdat_(char *str)
{
  time_t tim;
  tim = time(0);  
  strcpy(str,ctime(&tim));
}
