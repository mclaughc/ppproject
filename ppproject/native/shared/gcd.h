#pragma once
#include <cassert>

template <typename T>
static T gcd(T k, T m)
{
  assert(k >= T(0) && m >= T(0));
  while (k != m)
  {
    if (k > m)
      k -= m;
    else
      m -= k;
  }
  return k;
}
