

function tn

  j = SYSTIME(/JULIAN)
  CALDAT, j, M, D, Y, HH, MM, SS
  tn = string(format='(7(I4,I02,I02,"_",I02,I02,I02))', Y, M, D, HH, MM, SS)
  return, tn
end
