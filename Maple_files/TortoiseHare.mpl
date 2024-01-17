#
# The Tortoise and the Hare
# An algorithm for finding the period
# of an attracting cycle in an iteration
# a[n+1] = F(a[n]) of complex numbers
#
# Calling Sequence:
# (K,a,ratio) := TortoiseHare( a0, F, {tol=positive float}, {maxIt=posint} );
#
# Parameters: a[0] --- the initial point
#             F    --- the function F: C --> C
#             tol  --- the tolerance for deciding things are equal
#                      default 5.0e(-Digits) = Float(5,-Digits)
#             maxIt--- the maximum number of iterations to try, default 1000
#
# Method: as described in the Wikipedia article
#         https://en.wikipedia.org/wiki/Cycle_detection
#         with a modification for working over C and
#         using approximate equality, with a test for
#         contraction afterwards
#
# Output: k, period length K (positive integer)
#         a, first entry in the period (complex number, accurate to "tol")
#         r, ratio of closenesses on next pass (smaller than 1 is good)
#
# (c) Robert M. Corless 2023-09-10
# MIT License
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


TortoiseHare := proc( a0, F, {tol:=Float(5,-Digits)}, {maxIt:=1000} )
  local aHare, aTortoise, closeness, k;
  aTortoise := F(a0);
  aHare := F(aTortoise);
  if abs(aHare-aTortoise)<= tol then
    # Already converged
    return (1, aHare, abs(aHare-aTortoise)/abs(aTortoise-a0) );
  end if;
  for k to maxIt do
    aTortoise := F(aTortoise);
    aHare := F(aHare);
    aHare := F(aHare);
    closeness := abs( aHare - aTortoise );
    if closeness <= tol then
      break
    end if;
  end do;
  if k>maxIt then
    WARNING( "Maximum iterations exceeded." );
    return (maxIt,aHare,1.0);
  end if;
  # Post process: one more pass
  k := 0;
  aTortoise := aHare;
  for k to maxIt do
    aHare := F(aHare);
    if abs( aHare-aTortoise ) <= tol then
      return (k, aHare, abs(aHare-aTortoise)/closeness )
    end if;
  end do;
end proc:
