#ifndef SPINLOCK_H_INCLUDED
# define SPINLOCK_H_INCLUDED

typedef struct spinlock
{
  volatile unsigned long value;
} spinlock_t;

static inline void spinlock_init(spinlock_t* l)
{
  l->value = 0;
}

static inline int spinlock_trylock(spinlock_t* l)
{
  return __sync_bool_compare_and_swap(&l->value, 0, 1);
}

static inline void spinlock_lock(spinlock_t* l)
{
  while (!spinlock_trylock(l))
    ;
}

static inline void spinlock_unlock(spinlock_t* l)
{
  __sync_fetch_and_and(&l->value, 0);
}

#endif /* ! SPINLOCK_H_INCLUDED */
