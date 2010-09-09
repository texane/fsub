#ifndef CONFIG_H_INCLUDED
# define CONFIG_H_INCLUDED


#define CONFIG_KSIZE (20)
#define CONFIG_LSIZE (1 * CONFIG_KSIZE)
#define CONFIG_ASIZE (10 * CONFIG_LSIZE)
#define CONFIG_PAR_SIZE (0) /* 0 means 1 row per time */
#define CONFIG_ITER_COUNT 5
#define CONFIG_TIME 1
#define CONFIG_FSUB_SEQ 0
#define CONFIG_FSUB_GSL 1

/* pthread config */
#define CONFIG_FSUB_PTHREAD 0
#define CONFIG_THREAD_COUNT 8

/* xkaapi config */
#define CONFIG_FSUB_XKAAPI 1
#define CONFIG_PAR_GRAIN 1
#define CONFIG_SEQ_GRAIN 1


#endif /* ! CONFIG_H_INCLUDED */
