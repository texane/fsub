#ifndef ATOMIC_H_INCLUDED
# define ATOMIC_H_INCLUDED

static inline unsigned long read_atomic_ul(volatile unsigned long* ul)
{
  return *ul;
}

static inline void write_atomic_ul(volatile unsigned long* ul, unsigned long n)
{
  *ul = n;
}

static inline void inc_atomic_ul(volatile unsigned long* ul)
{
  __sync_fetch_and_add(ul, 1);
}

static inline void add_atomic_ul
(volatile unsigned long* ul, unsigned long n)
{
  __sync_fetch_and_add(ul, n);
}

static inline void dec_atomic_ul(volatile unsigned long* ul)
{
  __sync_fetch_and_add(ul, -1);
}

static inline int cas_atomic_ul
(volatile unsigned long* a, unsigned long b, unsigned long c)
{
  return __sync_bool_compare_and_swap(a, b, c);
}

#endif /* ! ATOMIC_H_INCLUDED */
