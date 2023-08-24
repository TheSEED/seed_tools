
use strict;
use Getopt::Long::Descriptive;
use LWP::UserAgent;

my @seed_list = ([PubSEED => 'http://pubseed.theseed.org/FIG'],
		 [PSEED => 'http://pseed.theseed.org'],
		 [CoreSEED => 'http://core.theseed.org/FIG'],
		 [OpenSEED => 'http://open.theseed.org/FIG'],
		 [RAST => 'http://rast.nmpdr.org'],
		);

my @names = map { $_->[0] } @seed_list;

my @options = (['output|o=s', 'Write output to the given file'],
	       ['pretty', 'Generate prettyprinted genome object'],
	       ['username=s' => 'RAST username'],
	       ['password=s' => 'RAST password'],
	       ['help|h', 'Print usage message and exit'],
	       [],
	       ["Valid SEED names are @names"],
	       );

my($opt, $usage) = describe_options("%c %o seed-name genome-id > output",
				    @options);

print($usage->text), exit if $opt->help;
print(STDERR $usage->text), exit 1 if @ARGV != 2;

my $seed = shift;
my $genome = shift;

my @url = grep { $_->[0] eq $seed } @seed_list;

if (@url != 1)
{
    print STDERR "Unknown SEED name $seed\n";
    print STDERR $usage->text;
    exit 1;
}

my $url = $url[0]->[1];

my $ua = LWP::UserAgent->new();

if ($seed eq 'RAST')
{
    $url .= "/genome_object.cgi?rast_job=$genome&username=" . $opt->username . "&password=" . $opt->password;
    $url .= "&pretty=1" if $opt->pretty;
}
else
{
    $url .= "/genome_object.cgi?genome=$genome";
    $url .= "&pretty=1" if $opt->pretty;
}

my $res = $ua->get($url);

if ($res->is_success)
{
    if ($opt->output)
    {
	open(F, ">", $opt->output) or die "Cannot write output " . $opt->output . ": $!\n";
	print F $res->content;
	close F;
    }
    else
    {
	print $res->content;
    }
}
else
{
    print STDERR "Error retrieving url $url: " . $res->content . "\n";
    exit 1;
}


		   
