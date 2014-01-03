function o(x) {
  print $1, $2, x
}

BEGIN {
  pft = 20
  lat = 2
  trees = 19
  grass = 18
  tr = 16
  te = 17
  bt = 15
  total = 14
  bne = 3
  bine = 4
  bns = 5
  tebs = 6
  ibs = 7
  tebe = 8
  trbe = 9
  tribe = 10
  trbr = 11
  bnee = 20
  trbee = 21
}

{
if ($trees > 2.5 && $trbe > .6*$trees) o(9)
else if ($trees > 2.5 && $trbr > .6*$trees) o(10)
else if ($trees > 2.5 && $tr > .5*$trees && (($trbe > $tebe && $trbe > $tebs) || ($trbr > $tebe && $trbr > $tebs))) o(8)
else if ($trees > 2.5 && $bt > .8*$trees && ($bnee > $bns || $ibs > $bns)) o(2)
else if ($trees > 2.5 && $bt > .8*$trees && ($bns > $bnee && $bns > $ibs)) o(1)
else if ($trees > 2.5 && $tebe > .5*$trees) o(6)
else if ($trees > 2.5 && $tebs > .5*$trees) o(5)
else if ($trees > 2.5 && $bt > .2*$trees) o(3)
else if ($trees > 2.5) o(7)
else if ($trees > .5 && $trees < 2.5 && $bt > .8*$trees && ($bnee > $bns || $ibs > $bns)) o(2)
else if ($trees > .5 && $trees < 2.5 && $bt > .8*$trees && ($bns > $bnee && $bns > $ibs)) o(1) 
else if ($trees > .5 && $trees < 2.5 && $trees > .8*$total) o(15)
else if ($trees > .5 && $trees < 2.5 && $total > 2.5) o(11)
else if ($trees > .5 && $trees < 2.5) o(12)
else if ($trees < .5 && $grass > .2 && $lat > 54) o(18)
else if ($grass > 2.0) o(13)
else if ($trees > .2 && $grass < 1.0) o(16)
else if ($grass > .2) o(14)
else if ($total > .2) o(16)
else if ($total <= .2) o(17)
}