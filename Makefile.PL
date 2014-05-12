use ExtUtils::MakeMaker;

use Config;
require 5.006;

$defines = $] < 5.008 ? "-DOLDPERL" : "-DNEWPERL";

my %options = %{
{
  NAME => 'Math::Complex_C',
  INC => '',
  DEFINE => $defines,
  LIBS => ['-lm'],
  LICENSE  => 'perl',
  VERSION_FROM => 'Complex_C.pm',
}
};
WriteMakefile(%options);

# Remove the Makefile dependency. Causes problems on a few systems.
sub MY::makefile { '' }
