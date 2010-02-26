#!/usr/bin/perl
# Count hits to the TINKER http site
# usage: webhits.pl files

# A hit is counted as one access to any TINKER page from a machine, max
# one per month. Multiple users on the same machine are lumped together,
# and same users on multiple machines are counted redundantly.

printf("\n");
printf("%s\n", ' Web Site Hits for the TINKER Molecular Modeling Package');
printf("\n");
use FileHandle;
$fh = new FileHandle;
$tothits = 0;
for $f (@ARGV)
  {
  my(%hit);
  if ($f =~ /\.gz$/) {$oc = "gzip -dcq $f|"}
  else {$oc = "<$f"};
  if ($fh->open($oc))
	{
	while(<$fh>)
	  {m%^(\S+).+GET /tinker/% && $hit{$1}++}
	$fh->close();
	}
  $hits = keys %hit;
  $tothits += $hits;
  printf("%8d %s\n", $hits, $f);
  }
printf("\n");
printf("%8d %s\n", $tothits, 'Total Hits');
