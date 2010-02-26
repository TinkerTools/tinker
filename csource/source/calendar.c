#include <time.h>

void
calendar(int* year,
         int* month,
         int* day, 
         int* hour,
         int* minute,
         int* second)
    {
    time_t now=time(NULL);
    struct tm* timestruct=localtime(&now);
    *year=timestruct->tm_year+1900;
    *month=timestruct->tm_mon+1;
    *day=timestruct->tm_mday;
    *hour=timestruct->tm_hour;
    *minute=timestruct->tm_min;
    *second=timestruct->tm_sec;
    }

int
calendar_(int* year,
          int* month,
          int* day, 
          int* hour,
          int* minute,
          int* second)
    {
    calendar(year,month,day,hour,minute,second);
    return 0;
    }
