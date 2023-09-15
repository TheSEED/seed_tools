#
# Copyright (c) 2003-2011 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
#
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License.
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#

package HTML;

use strict;
use Tracer;
use FIG;
use FIGRules;
use Carp;
use Data::Dumper;
use LWP::UserAgent;
use LWP::Simple;
use URI::Escape;  # uri_escape()
use URI::URL;
use HTTP::Request::Common;
use POSIX;
use CGI;

#use raelib; # now used for the excel function, that should eventually end up in here. Way too experimental!
my $raelib;


my $top_link_cache;


sub new
{
    my($class) = @_;

    my $self = {};

    return bless $self, $class;
}

sub top_link
{

    #
    # Determine if this is a toplevel cgi or one in one of the subdirs (currently
    # just /p2p).
    #

    return $top_link_cache if ($top_link_cache);

    my @parts = split(/\//, $ENV{SCRIPT_NAME});
    my $top;
    if (defined $parts[-2] && $parts[-2] eq 'FIG')
    {
        $top = '.';
#       warn "toplevel @parts\n";
    }
    elsif (defined $parts[-3] && $parts[-3] eq 'FIG')
    {
        $top = '..';
#       warn "subdir @parts\n";
    }
    else
    {
        $top = $FIG_Config::cgi_base;
#       warn "other @parts\n";
    }

    $top_link_cache = $top;
    return $top;
}

sub compute_html_header
{
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($additional_insert, $user, %options ) = @_;

    local $/ = "\n";

    my $header_name = $options{header_name} ? $options{header_name} : "html.hdr";
    my $tail_name = $options{tail_name} ? $options{tail_name} : "html.tail";

    my $html_hdr_file = "./Html/$header_name";
    if (! -f $html_hdr_file)
    {
        $html_hdr_file = "$FIG_Config::fig/CGI/Html/$header_name";
    }
    my @html_hdr = &FIG::file_read($html_hdr_file);

    # for my $k (sort keys %ENV) { warn "$k = $ENV{$k}\n"; }

    #
    # Determine if this is a toplevel cgi or one in one of the subdirs (currently
    # just /p2p).
    #

    my @parts = split(/\//, $ENV{SCRIPT_NAME});
    my $top;
    if ($parts[-2] eq 'FIG')
    {
        $top = '.';
#       warn "toplevel @parts\n";
    }
    elsif ($parts[-3] eq 'FIG')
    {
        $top = '..';
#       warn "subdir @parts\n";
    }
    else
    {
        $top = $FIG_Config::cgi_base;
#       warn "other @parts\n";
    }

    $options{no_fig_search} or push( @html_hdr, "<br><a href=\"$top/index.cgi?user=$user\">FIG search</a>\n" );

    if (@html_hdr)
    {
        my $insert_stuff;

        if (not $options{no_release_info})
        {
            my @ver = &FIG::file_head("$FIG_Config::fig_disk/CURRENT_RELEASE", 1);
            my $ver = $ver[0];
            chomp $ver;
            if ($ver =~ /^cvs\.(\d+)$/)
            {
                my $d = asctime(localtime($1));
                chomp($d);
                $ver .=  " ($d)";
            }
            my $host = &FIG::get_local_hostname();
            $insert_stuff = "SEED version <b>$ver</b> on $host";
        }

        if ($additional_insert)
        {
            $insert_stuff .= "<br>" . $additional_insert;
        }

        for $_ (@html_hdr)
        {
            s,(href|img\s+src)="/FIG/,$1="$top/,g;
                s,(\?user\=)\",$1$user",;
            if ($_ eq "<!-- HEADER_INSERT -->\n")
            {
                $_ = $insert_stuff;
            }
        }
    }

    return @html_hdr;
}

sub show_page {
    #warn "SHOWPAGE: cgi=", Dumper(@_);
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$html,$no_home, $alt_header, $css, $javasrc, $cookie, $options) = @_;
    my $i;
    Trace("Setting top link.") if T(3);
    my $top = top_link();

    # ARGUMENTS:
    #     $cgi is the CGI method
    #     $html is an array with all the html in it. It is just joined by "\n" (and not <br> or <p>
    #     $no_home eliminates ONLY the bottom FIG search link in a page
    #     $alt_header is a reference to an array for an alternate header banner that you can replace the standard one with
    #     $css is a reference to a hash. The key is the name of the CSS sheet and the value is the URL of that sheet. Note the usual rules about relative css urls
    #               the sheet named "Default" is considered to be the default style sheet, and if this is not set it points at $FIG_Config::HTML/css/default.css
    #               the sheet named "Sans Serif" is considered to the the first alternate, and if this is not set it points at $FIG_Config::HTML/css/sanserif.css
    #     $javasrc is a reference to an array of URLs to javascripts to be included (e.g. "FIG/Html/css/styleswitcher.js")
    #     $cookie is the name and value of the cookie to set. Note that you should probably use raelib->cookie to get/set your cookies
    #     $options is a reference to a hash of options that you can pass around the pages
    #
    # Find the HTML header
    #
    Trace("Reading tail.") if T(3);
    my $tail_name = $options->{tail_name} ? $options->{tail_name} : "html.tail";
    my $html_tail_file = "./Html/$tail_name";
    if (! -f $html_tail_file)
    {
        $html_tail_file = "$FIG_Config::fig/CGI/Html/$tail_name";
    }
    Trace("Extracting user name and header data.") if T(3);
    my $user = $cgi->param('user') || "";
    my @html_hdr;
    if ($alt_header && ref($alt_header) eq "ARRAY")
    {
       @html_hdr = @$alt_header;
    }
    else
    {
        @html_hdr = compute_html_header(undef,$user,%$options);
    }

    # RAE: I am offloading the handling of cookies to CGI.pm since I don't know how they are set up.
    # This modification adds the cookies if necessary

    # Note: 3/10/05 commented this line out pending the discussion of adding cookies into the seed that we are waiting to see about
    # to add cookies back in replace these two header lines with each other
    #my $hdr_thing = $cgi->header(-cookie=>$cookie);
    my $hdr_thing = $cgi->header();
    Trace("Printing HTML header: $hdr_thing.") if T(3);
    print $hdr_thing;
    Trace("Header printed.") if T(3);
    #
    #  The SEED header file goes immediately after <BODY>.  Figure out
    #  what parts of the HTML document skeleton are there, and fill in
    #  missing ones.
    #
    #  This list should be as comprehensive as feasible:
    #

    my %head_tag = ( base     => 1,
                     basefont => 1,
                     html     => 1,
                     isindex  => 1,
                     link     => 1,
                     meta     => 1,
                     nextid   => 1,
                     style    => 1,
                     title    => 1,
                   );

    #
    #  This list need not be comprehensive; it is just stopping conditions:
    #

    my %body_tag = ( a      => 1,
                     br     => 1,
                     center => 1,
                     form   => 1,
                     h1     => 1,
                     h2     => 1,
                     h3     => 1,
                     hr     => 1,
                     img    => 1,
                     p      => 1,
                     pre    => 1,
                     table  => 1
                   );

    my $html_line = -1;
    my $head_line = -1;
    my $meta_line = -1;
    my $base_line = -1;
    my $head_end_line = -1;
    my $body_line = -1;
    my $last_head_line = -1;  #  If no head tags are found, text goes at top.
    my $done = 0;
    Trace("Processing special cases.") if T(3);
    for ( $i = 0; $i < @$html; $i++ )
    {
        #  Some special cases:

        if ( $html->[$i] =~ /\<html\b/i )              { $html_line = $i }
        if ( $html->[$i] =~ /\<head\b/i )              { $head_line = $i }
        if ( $html->[$i] =~ /\<meta\s+html-equiv\b/i ) { $meta_line = $i }
        if ( $html->[$i] =~ /\<base\b/i )              { $base_line = $i }
        if ( $html->[$i] =~ /\<\/head\>/i )            { $head_end_line = $i }

        #  The content goes after this line:

        if ( $html->[$i] =~ /\<body[^0-9a-z]/i )
        {
            $body_line = $i;
            last;
        }

        #  Now the general case.
        #  Analyze all the html tags on the line:

        foreach ( $html->[$i] =~ /\<\/?([0-9a-z]+)/ig )
        {
            #  At first body tag, we stop the search and put the text
            #  after the last line with a head tag:

            if ( $body_tag{ lc $_ } )
            {
                $done = 1;
                last;
            }

            #  If this is a head tag, then move the marker forward

            elsif ( $head_tag{ lc $_ } )
            {
                $last_head_line = $i;
            }
        }
        last if $done;      # When done, break loop to avoid increment
    }

    #  Some sanity checks on structure:

    if ( 1 )
    {
        Trace("Sanity checks in progress.") if T(3);
        if ( $html_line >= 0 )
        {
            if ( ( $head_line >= 0 ) && ( $html_line > $head_line ) )
            {
                Trace("<HTML> tag follows <HEAD> tag.") if T(1);
            }
            if ( ( $head_end_line >= 0 ) && ( $html_line > $head_end_line ) )
            {
                Trace("<HTML> tag follows </HEAD> tag.") if T(1);
            }
        }
        if ( $head_line >= 0 )
        {
            if ( ( $head_end_line >= 0 ) && ( $head_line > $head_end_line ) )
            {
                Trace("<HEAD> tag follows </HEAD> tag.") if T(1);
            }
        }
    }
    Trace("Sanity checks complete.") if T(3);
    #
    #  Okay.  Let's put in the html header file, and missing tags:
    #
    #  <BODY> goes after last head line
    #
    #  RAE:
    #  Added the javascript for the buttons immediately after body.
    #  Note if no buttons are added we still (at the moment) add the script,
    #  but it only adds a little text (495 characters) to the html and noone will notice!
    #  RAE: This is now deprecated because everything is in an external file, FIG.js, included later
    if ( $body_line < 0 )
    {
        $body_line = $last_head_line + 1;
        Trace("Splicing body line at $body_line.") if T(3);
        splice( @$html, $body_line, 0, "<BODY>\n" );
    }

    #
    #  Seed page header (if it exists) goes after <BODY>
    #

    if (@html_hdr)
    {
        Trace("Splicing SEED page header after $body_line.") if T(3);
        splice( @$html, $body_line + 1, 0, @html_hdr );
    }

    #
    #  </HEAD> goes before <BODY>
    #

    if ( $head_end_line < 0 )
    {
        $head_end_line = $body_line;
        Trace("Splicing header terminater at $body_line.") if T(3);
        splice( @$html, $body_line, 0, "</HEAD>\n" );
    }

    # RAE:
    # Add css here
    # Note that at the moment I define these two sheets here. I think this should
    # be moved out, but I want to try it and see what happens.  css has the format:
    #
    # <link rel='stylesheet' title='default' href='/css/default.css' type='text/css'>
    Trace("Formatting CSS.") if T(3);
    # convert the default key to the right case. and eliminate dups
    foreach my $k (keys %$css) {if (lc($k) eq "default") {$css->{'Default'}=$css->{$k}}}

    if (!$css || !$css->{'Default'})
    {
       $css->{'Default'} = "Html/css/default.css";
    }
    if (!$css->{"Sans Serif"})
    {
       $css->{'Sans Serif'} = "Html/css/sanserif.css";
    }

    my $csstext = "<link rel='stylesheet' title='default' href='".$css->{'Default'}."' type='text/css'>\n";
    $csstext   .= "<link rel='alternate stylesheet' title='Sans Serif' href='".$css->{'Sans Serif'}."' type='text/css'>\n";

    foreach my $k (keys %$css)
    {
       next if (lc($k) eq "default" || lc($k) eq "sans serif");
       $csstext .= "<link rel='alternate stylesheet' title='$k' href='".$css->{$k}."' type='text/css'>\n";
    }

    $csstext   .= "<link rel='alternate'  title='SEED RSS feeds' href='Html/rss/SEED.rss' type='application/rss+xml'>\n";

    # RAE: also added support for external javascripts here.
    # we are cluttering the HTML code with all the javascripts when they could easily be in external files
    # this solution allows us to source other files

    # the file FIG.js contains most of the java script we use routinely. Every browser will just cache this and so
    # it will reduce our overhead.
    Trace("Formatting javascript.") if T(3);
    # $javasrc must be a ref to an array with urls (absolute or relative) to the javascripts
    push @$javasrc, "$FIG_Config::cgi_url/Html/css/FIG.js";
    foreach my $script (@$javasrc) {
        $csstext .= "<script src=\"$script\" type=\"text/javascript\"></script>\n";
    }

    Trace("Re-splicing the header terminator at $head_end_line.") if T(3);
    #  GJO 2011-09-13
    #  Changed to prefix to the </head> line, rather than replace it.
    $html->[$head_end_line] = $csstext . $html->[$head_end_line];

    #
    #  <BASE ...> goes before </HEAD>
    #

    if ( $base_line < 0 )
    {
        #
        #  Use a relative base address for pages.  Also, because I am
        #  worried about when FIG_config.pm gets updated (clean installs
        #  only, or every update?), I provide an alternative derivation
        #  from $cgi_url. -- GJO
        #
        # BASE href needs to be absolute. RDO.
        #
        #
#        $base_url = &FIG::cgi_url;
#       my $base_url = $FIG_Config::cgi_base;
#       if ( ! $base_url )                      # if cgi_base was not defined
#       {
#           $base_url = $FIG_Config::cgi_url;   # get the full cgi url
#           $base_url =~ s~^http://[^/]*~~;     # remove protocol and host
#           $base_url =~ m~/$~ || $base_url =~ s~$~/~; # check trailing slash
#       }

        $base_line = $head_end_line;
        #
        # RDO 2005-1006. Remove this so proxying works better.
        #
#        splice( @$html, $base_line, 0, "<BASE href=\"$base_url/\">\n" );
    }

    #
    #  <HTML> goes at the top of the output
    #

    if ( $html_line < 0 )
    {
        $html_line = 0;
        Trace("Splicing the HTML tag at $html_line.") if T(3);
        splice( @$html, $html_line, 0, "<HTML>\n" );
    }

    #
    #  <HEAD> goes after <HTML>
    #

    if ( $head_line < 0 )
    {
        $head_line = $html_line + 1;
        Trace("Splicing the HEAD tag at $head_line.") if T(3);
        splice( @$html, $head_line, 0, "<HEAD>\n" );
    }

    #
    #  <META html-equiv> goes after <HEAD>
    #

    if ( $meta_line < 0 )
    {
        $meta_line = $head_line + 1;
        Trace("Splicing the HEAD tag at $head_line.") if T(3);
        splice( @$html, $meta_line, 0, qq(<META http-equiv="Content-Type" content="text/html;charset=UTF-8" />\n) );
    }

    #
    #  Place FIG search link at bottom of page
    #
    Trace("Placing FIG search link.") if T(3);
    my @tail = -f $html_tail_file ? `cat $html_tail_file` : ();
    if (! $no_home)
    {
        my $user = $cgi->param('user') || "";
        push( @tail, "<hr><a href=\"index.cgi?user=$user\">FIG search</a>\n" );
    }

    #
    # See if we have a site-specific tail (for disclaimers, etc).
    #
    Trace("Placing site tail.") if T(3);
    my $site_tail = "$FIG_Config::fig_disk/config/site_tail.html";
    my $site_fh;
    if (open($site_fh, "<$site_tail"))
    {
        push(@tail, <$site_fh>);
        close($site_fh);
    }

    #
    #  Figure out where to insert The SEED tail.  Before </body>,
    #  or before </html>, or at end of page.
    #
    my @tags = ();
    Trace("Processing closing tags.") if T(3);
    for ($i=0; ($i < @$html) && ($html->[$i] !~ /\<\/body\>/i); $i++) {}
    if ($i >= @$html)        # </body> not found; look for </html>
    {
        push @tags, "\n</BODY>\n";
        # Even if tag is not found, index points to correct place for splice
        for ($i=0; ($i < @$html) && ($html->[$i] !~ /\<\/html\>/i); $i++) {}
        if ($i >= @$html)    # </html> not found; add it
        {
            push @tags, "</HTML>\n";
        }
    }

    if ( @tail )
    {
        Trace("Splicing tail.") if T(3);
        splice( @$html, $i, 0, @tail, @tags );
    }
    elsif ( @tags )
    {
        Trace("Splicing tags.") if T(3);
        splice( @$html, $i, 0, @tags );
    }

    Trace("Printing the HTML array.") if T(3);
    # RAE the chomp will return any new lines at the ends of elements in the array,
    # and then we can join  with a "\n". This is because somethings put newlines in,
    # and others don't. This should make nicer looking html
    #
    # chomp(@$html);
    # print join "\n", @$html;
    #
    # Apparently the above still breaks things. This is the correct code:
    # Removed superfluous assignment to variable: GJO -- 2011-08-13
    #
    foreach ( @$html )
    {
        if (T(4)) {
            my $escapedLine = Tracer::Clean($_);
            Trace("Printing:\n$escapedLine") if T(4);
        }
        print $_;
    }

}


=head1 make_table

The main method to convert an array into a table.

The col_hdrs are set to the <th> headers, the $tab is an array of arrays. The first is the rows, and the second is the columns. The title is the title of the table.

The options define the settings for the table such as border, width, and class for css formatting. If the option "excelfile" is set to a filename that will be created in FIG_Config::temp, and the link included that allows the user to download the file as an excel file.

=cut

sub make_table {
    my($col_hdrs,$tab,$title, %options ) = @_;
    my(@tab);

    my $border = defined $options{border} ? "border=\"$options{border}\"" : "border";
    my $width = defined $options{width} ? "width=\"$options{width}\"" : "";
    my $class = defined $options{class} ? "class=\"$options{class}\"" : "";
    push( @tab, "\n<table $border $width $class>\n",
                "\t<caption><b>$title</b></caption>\n",
                "\t<tr>\n\t\t"
              . join( "\n", map { &expand($_, "th") } @$col_hdrs )
              . "\n\t</tr>\n"
        );
    my($i);

    my $row;
    foreach $row (@$tab)
    {
        push( @tab, "\t<tr>\n"
                  . join( "\n", map { &expand($_) } @$row )
                  . "\n\t</tr>\n"
            );
    }
    push(@tab,"</table>\n");

    # excelfile should be appropriate for a filename (no spaces/special characters)
    if (defined $options{"excelfile"}) {
        if (! defined($raelib)) {
            require raelib;
            $raelib = new raelib;
        }
        push @tab, $raelib->tab2excel($col_hdrs,$tab,$title,\%options,$options{"excelfile"})}

    return join("",@tab);
}

sub abstract_coupling_table {
    my($cgi,$prot,$coupling) = @_;
    my %fc;

    my $col_hdrs = ["coupled to","Score","Type of Coupling", "Type-specific Data"];
    my $tab = [];
    my %by_peg;
    foreach my $x (@$coupling)
    {
        my($peg2,$psc,$type,$extra) = @$x;
        if (($type !~ /^[ID]FC$/) || (! $fc{$peg2}))
        {
            if ($type =~  /^[ID]FC$/)
            {
                $fc{$peg2} = 1;
            }

            $by_peg{$peg2} += $psc;
        }
    }

    foreach my $x (sort { ($by_peg{$b->[0]} <=> $by_peg{$a->[0]})
                          or ($a->[0] cmp $b->[0])
                          or ($b->[1] <=> $a->[1])
                          or ($a->[2] cmp $b->[2]) } @$coupling)
    {
        my($peg2,$psc,$type,$extra) = @$x;
        push(@$tab,[&fid_link($cgi,$peg2,1),$psc,$type,&set_prot_links($cgi,join(", ",@$extra))]);
    }


     my $help = "<a href=\"Html/abstract_coupling.html\" target=\"SEED_or_SPROUT_help\">for help</a>";
#    my @html = &make_table($col_hdrs,$tab,"Abstract Coupling Data for $prot");
#    push(@html,"<hr>\n",$cgi->h3($help),"<br>");
#    return @html;

    return &make_table($col_hdrs,$tab,"Abstract Coupling Data for $prot [$help]");
}

sub expand {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my( $x, $tag ) = @_;

    $tag = "td" unless $tag;
    my $endtag = $tag;

    # RAE modified this so that you can pass in a reference to an array where
    # the first element is the data to display and the second element is optional
    # things like colspan and align. Note that in this case you need to include the td
    # use something like ["some data to appear", "td colspan=4 bgcolor=gray"]

    # per GJO's request modified this line so it can take any tag.
    if ( ref($x) eq "ARRAY" ) { ($x, $tag) = @$x; $tag =~ /^(\S+)/; $endtag = $1 }

    if ( $x =~ /^\@([^:]+)\:(.*)$/ )
    {
        return "\t\t<$tag $1>$2</$endtag>";
    }
    else
    {
        return "\t\t<$tag>$x</$endtag>";
    }
}


=head2 merge_table_rows()

Merge table rows together. This will merge a table so that adjacent cells with the same content will only be shown once.

Something like this:

    -----------------------
    |    1     |    a     |
    -----------------------
    |    1     |    b     |
    -----------------------
    |    2     |    c     |
    -----------------------
    |    3     |    d     |
    -----------------------
    |    4     |    d     |
    -----------------------
    |    5     |    d     |
    -----------------------

Will become:

    -----------------------
    |          |    a     |
    |    1     |-----------
    |          |    b     |
    -----------------------
    |    2     |    c     |
    -----------------------
    |    3     |          |
    ------------          |
    |    4     |    5     |
    ------------          |
    |    5     |          |
    -----------------------


The method takes two arguments. The reference to the array that is the table ($tab). This is the standard table that is created for HTML.pm to draw, and a reference to a hash of columns that you don't want to merge together. The reference to the hash is optional, and if not included, everything will be merged.

 $tab=&HTML::merge_table_rows($tab);

 or

 $skip=(1=>1, 3=>1, 5=>1);
 $tab=&HTML::merge_table_rows($tab, $skip);  # will merge all columns except 1, 3 and 5. Note the first column in the table is #0


=cut




sub merge_table_rows {
 # RAE:
 # Experimental piece of code. We often want to have rows or columns were cells are merged. It just looks so much nicer
 # this block should merge adjacent rows that have the same text in them.
 # use like this:
 #      $tab=&HTML::merge_table_rows($tab);
 # before you do a make_table call

 my $self=shift if UNIVERSAL::isa($_[0],__PACKAGE__);
 my ($tab, $skip)=@_;

 my $newtable;
 my $lastrow;
 my $rowspan;
 my $refs;

 for (my $y=0; $y <= $#$tab; $y++) {
 #$y is the row in the table;
  for (my $x=0; $x <= $#{$tab->[$y]}; $x++) {
   # this is the user definable columns not to merge
   if ($skip->{$x})
   {
    $newtable->[$y]->[$x] = $tab->[$y]->[$x];
    next;
   }

   #$x is the column in the table
   # if the column in the row we are looking at is the same as the column in the previous row, we don't add
   # this cell to $newtable. Instead we increment the rowspan of the previous row by one

   # handle cells that are references to arrays
   if (ref($tab->[$y]->[$x]) eq "ARRAY") {$refs->[$y]->[$x]=$tab->[$y]->[$x]->[1]; $tab->[$y]->[$x]=$tab->[$y]->[$x]->[0]}

   # now we go back through the table looking where to draw the merge line:
   my $lasty=$y;
   while ($lasty >= 0 && $tab->[$y]->[$x] eq $tab->[$lasty]->[$x]) {$lasty--}
   $lasty++; # this is the last identical cell. If lasty==y it is the current cell, so we just save the data. Otherwise we increment the rowspan
   if ($lasty == $y) {
    # we always want to have something in rows that may otherwise be empty but should be there (see below)
    unless ($tab->[$y]->[$x]) {$tab->[$y]->[$x]=" &nbsp; "}
    $newtable->[$y]->[$x] = $tab->[$y]->[$x];
   }
   else {$rowspan->[$lasty]->[$x]++}
  }
 }

 # now just join everything back together
 for (my $y=0; $y <= $#$tab; $y++) {
  for (my $x=0; $x <= $#{$tab->[$y]}; $x++) {
   if ($rowspan->[$y]->[$x]) {
    if ($refs->[$y]->[$x]) {$refs->[$y]->[$x] .= " rowspan=". ($rowspan->[$y]->[$x]+1)}
    else {$refs->[$y]->[$x] = "td rowspan=". ($rowspan->[$y]->[$x]+1)}
    $newtable->[$y]->[$x]=[$newtable->[$y]->[$x], $refs->[$y]->[$x]];
   }
   elsif ($newtable->[$y]->[$x] && $refs->[$y]->[$x]) {
    $newtable->[$y]->[$x]=[$newtable->[$y]->[$x], $refs->[$y]->[$x]];
   }
  }
 }


 # finally we have to remove any completely empty cells that have been added by the array mechanism
 # (e.g. if you define $a->[2] then $a->[0] and $a->[1] are now undef).
 # that is why in the loop above I replace empty cells with nbsp. They are now not undef!
 # I am sure that Gary can do this in one line, but I am hacking.
 my @trimmed;
 foreach my $a (@$newtable) {
  my @row;
  foreach my $b (@$a) {
   push @row, $b if ($b);
  }
  push @trimmed, \@row;
 }

 return \@trimmed;
}




sub set_ec_links {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$x) = @_;
    my($before,$match,$after);

    if ($x =~ /^(.*)(EC \d+\.\d+\.\d+\.\d+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_ec_links($cgi,$before) . &HTML::ec_link($match) . &set_ec_links($cgi,$after);
    }
    return $x;
}

sub ec_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($role) = @_;

    if ($role =~ /(\d+\.\d+\.\d+\.\d+)/)
    {
        return "<a href=\"http://www.genome.ad.jp/dbget-bin/www_bget?ec:$1\">$role</a>";
    }
    else
    {
        return $role;
    }
}

sub role_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$role) = @_;

    my $roleR = ($role =~ /^(\d+\.\d+\.\d+\.\d+)\s+-\s+/) ? $1 : $role;
    my $user = $cgi->param('user');
    if (! $user) { $user = "" }
    my $link = $cgi->url() . "?role=$roleR&user=$user";
    $link =~ s/[a-z]+\.cgi\?/pom.cgi?/;
    return "<a href=$link>$role</a>";
}

=head2 fid_link

Get a link to a fid.

use: my $html=&HTML::fid_link($cgi, $fid, Local, Just_URL, Full_Path);

Local is a boolean means to eliminate the fig|genome_id.peg. from the text of the link.

Just_URL will only return the URL and not the HTML code. The default is to return the full code.

Full_Path is a boolean that will get the full path to the URL not just a relative path. This is required in pages where the base href changes (e.g. if an image is imported like on the metabolic pages).

=cut


sub fid_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$fid,$local,$just_url,$fullpath) = @_;
    Trace("Creating link for feature $fid.") if T(4);
    my $err=join(" ", $cgi,$fid,$local,$just_url,$fullpath);

    my @params = ( '' );
    push @params, 'user=' . $cgi->param('user') if $cgi->param('user');

    my $text;
    if ( ! $just_url ) {
        $text = $fid;
        # If we're in local mode, we remove everything but the final number from the fig ID.
        ( $text =~ s/^.*\.peg\.(\d+)$/$1/ ) || ( $text =~ s/^.*\.([^.]+\.\d+)$/$1/ ) if $local;
    }

    if ( $FIG_Config::seedviewer_enabled ) {
        # We're NMPDR. Compute the link to the seed viewer feature page.
        my $link = "seedviewer.cgi?page=Annotation&feature=$fid" . join( '&', @params );

        if ($fullpath) {
            # Full-path mode: add the base URL.
            $link = "$FIG_Config::cgi_url/$link";
        }

        return $just_url ? $link : "<a href=\"$link\">$text</a>";
    }

    my ( $type, $num ) = $fid =~ /^fig\|\d+\.\d+\.([a-zA-Z]+)\.(\d+)/;
    if ( $type && $num )
    {
        my $top = $fullpath ? $FIG_Config::cgi_url : top_link();

        my $link;
        my $sprout = $cgi->param('SPROUT') || 0;
        push @params, 'SPROUT=1'                            if $sprout;
        push @params, '48hr_job=' . $cgi->param("48hr_job") if $cgi->param("48hr_job");
        push @params, 'new_framework=1'                     if $cgi->param('new_framework');
        if ($type ne "peg" && ! $sprout)
        {
           Trace("Creating feature link for $fid.") if T(4);
           $link = "$top/feature.cgi?feature=$fid" . join( '&', @params );
        }
        else
        {
            Trace("Creating protein link for $fid.") if T(4);
            push @params, 'translate=1' if $cgi->param('translate');

### This used to be
###     my $link = &FIG::cgi_url . "/protein.cgi?prot=$fid&user=$user$trans$sprout";
###
### The cost became prohibitive in the subsystem spreadsheets.  Hence, we cache the value
###
### RAO

            #if (! $cgi_url) { $cgi_url = &FIG::cgi_url }
            #$link = $cgi_url . "/protein.cgi?prot=$fid&user=$user$trans$sprout";
            $link = "$top/protein.cgi?prot=$fid" . join( '&', @params );
        }

        return $just_url ? $link : "<a target=_blank href='$link'>$text</a>";
    }
    return $fid;
}

sub family_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($family,$user) = @_;

    return $family;
}

=head2 evidence_codes_explain

Given an evidence code, returns a string that explains this eveidence code.

=cut

sub evidence_codes_explain {
 my($ec)=@_;
 return unless ($ec);

 $ec=uc($ec);
 return "IDA: Inferred from Direct Assay" if ($ec =~ /^IDA/);
 return "IGI: Inferred from Genetic Interaction" if ($ec =~ /^IGI/);
 return "TAS: Traceable Author Statement" if ($ec =~ /^TAS/);
 return "ISU: in subsystem unique" if ($ec =~ /^ISU/);
 return "$ec: in subsystem duplicates" if ($ec =~ /^IDU/);
 return "$ec: in cluster with" if ($ec =~ /^ICW/);
 return "FF: in FIGfam" if ($ec =~ /^FF/);
 return "CWN: clustered with nonhypothetical" if ($ec =~ /^CWN/);
 return "CWH: clustered, but only with hypotheticals" if ($ec =~ /^CWH/);
 return "DLIT: literature references to this gene exist" if ($ec =~ /^DLIT/);
 return "ILIT: no references to this gene exist, but they do to other genes with the same functional role" if ($ec =~ /^ILIT/);
 return "$ec: unknown!";
}

sub get_html {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my( $url, $type, $kv_pairs) = @_;
    my( $encoded, $ua, $args, @args, $out, @output, $x );

    $ua = new LWP::UserAgent;
    $ua->timeout( 900 );
    if ($type =~/post/i)
    {
        $args = [];
        foreach $x (@$kv_pairs)
        {
            push(@$args, ( $x->[0], $x->[1]) );
        }
        my $request  = POST $url, $args;
        my $response = $ua->request($request);
        $out = $response->content;
    }

    if ($type =~/get/i)
    {
        @args = ();
        foreach $x (@$kv_pairs)
        {
            push( @args, "$x->[0]=" . uri_escape($x->[1]) );
        }

        if (@args > 0)
        {
            $url .= "?" . join("&",@args);
        }
        my $request = new HTTP::Request('GET', $url);
        my $response = $ua->request($request);

        if ($response->is_success)
        {
            $out = $response->content;
        }
        else
        {
            $out = "<H1>Error: " . $response->code . "</H1>" . $response->message;
        }
    }
#   set up a document with proper eol characters
    @output = split(/[\012\015]+/,$out);
    foreach $out (@output) { $out .= "\n"; }

#   Now splice in a line of the form <base href=URL> to cause all relative links to work
#   properly.  Remove the header.
    my $i;
    for ($i=0; ($i < @output) && ($output[$i] !~ /^\s*\</); $i++) {}
    if ($i < @output) {
        splice(@output,0,$i);
    }

    for ($i=0; ($i < @output) && ($output[$i] !~ /\<body\>/i); $i++) {}
    if ($i == @output)
    {
        $i = -1;
    }
    splice(@output,$i+1,0,"<base href=\"$url\">\n");
    return @output;
}

sub trim_output {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($out) = @_;
    my ($i, $j);

    for ($i=0; ($i < @$out) && ($out->[$i] !~ /^\</); $i++) {}
    splice(@$out,0,$i);

    for ($i=0; ($i < @$out) && ($out->[$i] !~ /\<body\>/i); $i++) {}
    if ($i == @$out)
    {
        for ($i=0; ($i < @$out) && ($out->[$i] !~ /\<html\>/i); $i++) {}
        if ($i == @$out)
        {
            $i = -1;
        }
    }
    for ($j=$i+1; ($j < @$out) && ($out->[$j] !~ /^\<hr\>$/); $j++) {}
    if ($j < @$out)
    {
        splice(@$out,$i+1,($j-$i));
    }

    for ($i=0; ($i < @$out) && ($out->[$i] !~ /\<\/body\>/i); $i++) {}
    if ($i == @$out)
    {
        for ($i=0; ($i < @$out) && ($out->[$i] !~ /\<\/html\>/i); $i++) {}
    }

    for ($j=$i-1; ($j > 0) && ($out->[$j] !~ /FIG search/); $j--) {}
    if ($j > 0)
    {
        #
        # Hm. We would have tried using the options here:
        # my $tail_name = $options{tail_name} ? $options{tail_name} : "html.tail";
        # but they're not passed in. So use the default html.tail.
        #
        my $html_tail_file = "./Html/html.tail";
        my @tmp = `cat $html_tail_file`;
        my $n = @tmp;
        splice(@$out,$j-$n,$n+1);
    }
}

=head2 alias_rl

Returns the url that links to an external page showing information about the given alias.
The type of the alias will be determined by the prefix (i.e. 'tr|' for Trembl) If the type
cannot be determined, the function will return undef.

use: my $html=&HTML::alias_url($alias, $type);

=cut

sub alias_url {
  shift if UNIVERSAL::isa($_[0],__PACKAGE__);

  my ($id, $type) = @_;
  
  if ($type eq "SEED") { # 1
    return "http://seed-viewer.theseed.org/linkin.cgi?id=$id";
  }
  elsif ($type eq "UniProt") {
    return "http://www.uniprot.org/entry/$id";
  }
  elsif ($type eq "UniProt_ac") { # 2
    return "http://www.uniprot.org/entry/$id";
  }
  elsif ($type eq "UniProt_id") { # 3
    return "http://www.uniprot.org/entry/$id";
  }
  elsif ($type eq "EntrezGene") { # 4
    return "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id=$id";
  }
  elsif ($type eq "RefSeq") { # 5
    return "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id=$id";
  }
  elsif ($type eq "GIID") { # 6
    return "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id=$id";
  }
  elsif ($type eq "NCBI") {
    return "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id=$id";
  }
  elsif ($type eq "PDB") { # 7
    $id =~ s/\:\w//;
    return "http://www.rcsb.org/pdb/explore/explore.do?structureId=$id";
  }
  elsif ($type eq "PFAM") { # 8
    return "http://pfam.janelia.org/family?acc=$id";
  }
  elsif ($type eq "GO") { # 9
    return "http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=$id";
  }
  elsif ($type eq "PIRSF") { # 10
    return "http://pir.georgetown.edu/cgi-bin/ipcSF?id=$id";
  }
  elsif ($type eq "IPI") { # 11
    return "http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-newId+[IPI-AllText:".$id."*]+-lv+30+-view+SeqSimpleView+-page+qResult";
  }
  elsif ($type eq "UniRef_100") { # 12
    return "http://www.uniprot.org/entry/$id";
  }
  elsif ($type eq "UniRef_90") { # 13
    return "http://www.uniprot.org/entry/$id";
  }
  elsif ($type eq "UniRef_50") { # 14
    return "http://www.uniprot.org/entry/$id";
  }
  elsif ($type eq "UniParc") { # 15
    return "http://www.uniprot.org/entry/$id";
  }
  elsif ($type eq "PIR-PSD") { # 16
    return "http://pir.georgetown.edu/cgi-bin/pir_psd_get.pl?id=$id";
  }
  elsif ($type eq "Taxon_ID") { # 17
    return "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=$id";
  }
  elsif ($type eq "OMIM") { # 18
    return "http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=$id";
  }
  elsif ($type eq "UniGene") { # 19
    return "http://www.ncbi.nlm.nih.gov/sites/entrez?db=unigene&cmd=search&term=$id";
  }
  elsif ($type eq "Ensemble_ID") { # 20
    #return "$id";
  }
  elsif ($type eq "PMID") { # 21
    return "http://www.ncbi.nlm.nih.gov/pubmed/$id";
  }
  elsif ($type eq "EMBL_DNA_AC") { # 22
    return "http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-e+[EMBL:".$id."]+-newId";
  }
  elsif ($type eq "EMBL_Protein_AC") { # 23
    $id =~ s/\.\d//;
    return "http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-e+[{EMBL}-ProteinID:".$id."]";
  }
  elsif ($type eq "CMR") { # 24
    if ($id =~ /^\d+$/) {
      return "http://cmr.jcvi.org/cgi-bin/CMR/shared/GenePage.cgi?type=PID&acc=".$id;
    } else {
      return "http://cmr.jcvi.org/tigr-scripts/CMR/shared/GenePage.cgi?locus=".$id;
    }
  }
  elsif ($type eq "SwissProt"){ # 25
      return "http://www.uniprot.org/entry/$id";
  }
  elsif ($type eq "IMG"){ # 26
      return "http://img.jgi.doe.gov/cgi-bin/pub/main.cgi?page=geneDetail&gene_oid=$id";
  }
  elsif ($type eq "KEGG") { # 27
      my ($pre,$post) = ($id) =~ /(.*):(.*)/;
      return "http://www.genome.ad.jp/dbget-bin/www_bget?" . $pre . "+" . $post;
  }
  elsif ($type eq "ASAP") { # 28
    return "https://asap.ahabs.wisc.edu/asap/feature_info.php?FeatureID=$id";
  }
  elsif ($type eq "GenBank") { # 29
      return "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id=$id";
  }
  elsif ($type eq "DBJ") { # 30
      return "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id=$id";
  }
  elsif ($type eq "SCOP") { # 31
      $id =~ s/\:\w//;
      return "http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?pdb=$id";
  }
  elsif ($type eq "CATH") { # 32
      $id =~ s/\:\w//;
      return "http://www.cathdb.info/cgi-bin/CATHSrch.pl?type=PDB&query=$id";
  }
  elsif ($type eq "FSSP") { # 33
      $id =~ s/\:\w//;
      return "http://ekhidna.biocenter.helsinki.fi/dali/daliquery?find=$id";
  }
  elsif ($type eq "MMDB") { # 34
      $id =~ s/\:\w//;
      return "http://www.ncbi.nlm.nih.gov/Structure/mmdb/mmdbsrv.cgi?form=6&db=tDopt=s&uid=$id";
  }
  elsif ($type eq "PDBsum") { # 35
      $id =~ s/\:\w//;
      return "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode=$id";
  }

  return undef;
}

sub set_prot_links {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$x) = @_;
    my($before,$match,$after);

    if ($x =~ /^(.*)(fig\|\d+\.\d+\.peg\.\d+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::fid_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)\b([NXYZA][PM]_[0-9\.]+)\b(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::refseq_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(gi\|\d+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::gi_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(gb\|\S+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::gb_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(img\|\d+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::img_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(tigr\|\w+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::tigr_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(tigrcmr\|\w+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::tigr_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(cmr\|\w+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::cmr_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(dbj\|\S+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::dbj_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)\b(eric\|\S+)\b(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::eric_link($cgi,$match) . &set_prot_links($cgi,$after);
    }

    elsif ($x =~ /^(.*)\b(bhb\|.*?)\b(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::bhb_link($cgi,$match) . &set_prot_links($cgi,$after);
    }

    elsif ($x =~ /^(.*)\b(apidb\|[0-9\.a-z_]+)\b(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::apidb_link($cgi,$match) . &set_prot_links($cgi,$after);
    }

    elsif ($x =~ /^(.*)\b(patric\|.*?)\b(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::patric_link($cgi,$match) . &set_prot_links($cgi,$after);
    }

    elsif ($x =~ /^(.*)\b(vbrc\|.*?)\b(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::vbrc_link($cgi,$match) . &set_prot_links($cgi,$after);
    }

    elsif ($x =~ /^(.*)\b(vectorbase\|.*?)\b(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::vectorbase_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~  /^(.*)(uni\|[A-Z0-9]{6})(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::uni_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(sp\|[A-Z0-9]{6})(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::sp_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(pirnr\|NF\d+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::pir_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(kegg\|[a-z]{2,4}:[a-zA-Z_0-9]+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::kegg_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(Ensembl[a-zA-Z]+:[a-zA-Z_0-9\.]+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::ensembl_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(EntrezGene:[a-zA-Z_0-9\.]+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::entrezgene_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(MIM:[a-zA-Z_0-9\.]+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::mim_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(HGNC:[a-zA-Z_0-9\.]+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::hgnc_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(UniGene:[a-zA-Z_0-9\.]+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::unigene_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
# IPI stopped working. turn off for now.
#    elsif ($x =~ /^(.*)(IPI:[a-zA-Z_0-9\.]+)(.*)/s)
#    {
#        $before = $1;
#        $match = $2;
#        $after = $3;
#        return &set_prot_links($cgi,$before) . &HTML::ipi_link($cgi,$match) . &set_prot_links($cgi,$after);
#    }
    elsif ($x =~ /^(.*)(WP:[a-zA-Z_0-9\.]+)(.*)/s)
    {
        #wormbase

        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::wp_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(FB:[a-zA-Z_0-9\.]+)(.*)/s)
    {
        #flybase

        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::fb_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(FlyBaseORFNames:[a-zA-Z_0-9\.]+)(.*)/s)
    {
        #flybase

        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::fborf_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(SGD_LOCUS:[a-zA-Z_0-9\.]+)(.*)/s)
    {
        #flybase

        $before = $1;
        $match = $2;
        $after = $3;
        return &set_prot_links($cgi,$before) . &HTML::sgd_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    elsif ($x =~ /^(.*)(tr\|[a-zA-Z0-9]+)(.*)/s)
    {

      $before = $1;
      $match = $2;
      $after = $3;

      return &set_prot_links($cgi,$before) .  &HTML::trembl_link($cgi,$match) . &set_prot_links($cgi,$after);
    }
    return $x;
}

sub lit_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($id) = @_;

    if ( $id =~ /([di]lit)\((\d+)\)/ )
    {
	my($pre, $pubmed_id) = ($1, $2);
	my $pubmed_url = &alias_url($pubmed_id, 'PMID');
	$id = "$pre(<a href='$pubmed_url' target=_blank>$pubmed_id</a>)";
    }
    return $id;
}

sub refseq_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$id) = @_;

    if ($id =~ /^[NXYZA]P_/)
    {
        return "<a href='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=protein&cmd=search&term=$id' target=_blank>$id</a>";
    }
    elsif ($id =~ /^[NXYZA]M_/)
    {
        return "<a href='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=nuccore&cmd=search&term=$id' target=_blank>$id</a>";
    }
}

sub gi_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$gi) = @_;

    if ($gi =~ /^gi\|(\d+)$/)
    {
        return "<a href='http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=Protein&list_uids=$1&dopt=GenPept' target=_blank>$gi</a>";
    }
    return $gi;
}

sub gb_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$gb) = @_;

    if ($gb =~ /^gb\|(\S+)$/)
    {
        return "<a href='http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id=$1' target=_blank>$gb</a>";
    }
    return $gb;
}

sub tigr_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$tigr) = @_;

    if ($tigr =~ /^tigr\|(NT|ntbp|ntbpA|BA|BMAA|BXB|GBA)(\w+)$/)
    {
        my $id=$1.$2;
        return "<a href=\"http://pathema.tigr.org/tigr-scripts/pathema/shared/GenePage.cgi?locus=$id\" target=_blank>$tigr</a> (Pathema)";
    }
    elsif ($tigr =~ /^tigr(cmr)?\|(\S+)$/)
    {
        return "<a href=\"http://www.tigr.org/tigr-scripts/CMR2/GenePage.spl?locus=$2\" target=_blank>$tigr</a>";
    }
    return $tigr;
}

sub cmr_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$cmr) = @_;

    if ($cmr =~ /^cmr\|(\d+)$/)
    {
        my $id=$1;
	return "<a href=\"http://cmr.jcvi.org/cgi-bin/CMR/shared/GenePage.cgi?type=PID&acc=".$id."\" target=_blank>$cmr</a>";
    }
    elsif ($cmr =~ /^cmr\|(\S+)$/)
    {
	my $id = $1;
	return "<a href=\"http://cmr.jcvi.org/tigr-scripts/CMR/shared/GenePage.cgi?locus=".$id."\" target=_blank>$cmr</a>";
    }
    return $cmr;
}

sub eric_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$eric) = @_;

    if ($eric =~ /^eric\|(\S+)/)
    {
        return "<a href=\"https://asap.ahabs.wisc.edu/asap/feature_info.php?FeatureID=$1\" target=_blank>$eric</a>";
    }
    return $eric;
}

sub dbj_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$dbj) = @_;

    if ($dbj =~ /^dbj\|(\S+)/)
    {
        return "<a href=\"http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id=$1\" target=_blank>$dbj</a>";
    }
    return $dbj;
}

sub bhb_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$bhb) = @_;

    return "<a href=\"http://www.biohealthbase.org\" target=_blank>$bhb</a>";
}

sub apidb_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$api) = @_;

    if ($api =~ /apidb\|(.*?)\.(.*)$/)
    {
        return "<a href=\"http://www.apidb.org/cgi-bin/redirect.cgi?taxon_id=$1&source_id=$2\" target=_blank>$api</a>";
    }
    return $api;
}

sub patric_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$patric) = @_;

    if ($patric =~ /patric\|(.*)/)
    {
        return "<a href=\"https://patric.vbi.vt.edu/software/curationTool/gep/pgiCuration.php?locus_name=$1\" target=_blank>$patric</a>";
    }
    return $patric;
}

sub vbrc_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$vbrc) = @_;

    if ($vbrc =~ /vbrc\|(.*)/)
    {
        return "<a href=\"http://www.biovirus.org/gene_detail.asp?name=$1\" target=_blank>$vbrc</a>";
    }
    return $vbrc;
}

sub vectorbase_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$vec) = @_;
    return "<a href=\"http://www.vectorbase.org\" target=_blank>$vec</a>";
}

sub uniprot_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$uni) = @_;

    if ($uni =~ /^(sp|tr|uni)\|(\S+)$/)
    {
	return "<a href='" . &HTML::alias_url($2, 'UniProt') . "' target=_blank>$uni</a>";
    }
    return $uni;
}

sub uni_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$uni) = @_;

    if ($uni =~ /^uni\|(\S+)$/)
    {
        #return "<a href=http://www.pir.uniprot.org/cgi-bin/upEntry?id=$1>$uni</a>";
        #return "<a href='http://www.ebi.uniprot.org/uniprot-srv/uniProtView.do?proteinAc=$1' target=_blank>$uni</a>";
	return "<a href='" . &HTML::alias_url($1, 'UniProt') . "' target=_blank>$uni</a>";
    }
    return $uni;
}

sub sp_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$sp) = @_;

    if ($sp =~ /^sp\|(\S+)$/)
    {
        #return "<a href='http://us.expasy.org/cgi-bin/get-sprot-entry?$1' target=_blank>$sp</a>";
	return "<a href='" . &HTML::alias_url($1, 'UniProt') . "' target=_blank>$sp</a>";
    }
    return $sp;
}

sub trembl_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$tr) = @_;

    if ($tr =~ /^tr\|(\S+)$/)
    {
        #return "<a href='http://us.expasy.org/cgi-bin/get-sprot-entry?$1' target=_blank>$tr</a>";
	return "<a href='" . &HTML::alias_url($1, 'UniProt') . "' target=_blank>$tr</a>";
    }
    return $tr;
}

sub pir_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$pir) = @_;

    if ($pir =~ /^pirnr\|(NF\d+)$/)
    {
        return "<a href='http://pir.georgetown.edu/cgi-bin/nfEntry.pl?id=$1' target=_blank>$pir</a>";
    }
    return $pir;
}

sub kegg_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$kegg) = @_;

    if ($kegg =~ /^kegg\|([^:]+):(\S+)$/)
    {
        return "<a href='http://www.genome.ad.jp/dbget-bin/www_bget?$1+$2' target=_blank>$kegg</a>";
    }
    return $kegg;
}

sub img_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$img) = @_;

    if ($img =~ /^img\|(\S+)$/)
    {
        return "<a href='http://img.jgi.doe.gov/cgi-bin/pub/main.cgi?page=geneDetail&gene_oid=$1' target=_blank>$img</a>";
    }
    return $img;
}

sub ensembl_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$ensembl) = @_;

    if ($ensembl =~ /^(\S+):(\S+)$/)
    {
        my $what=$1;
        my $key=$2;
        my $idx="All";
        if ($what eq "EnsemblGene") { $idx = "Gene" }
        if ($what eq "EnsemblTranscript") { $idx = "All" }
        if ($what eq "EnsemblProtein") { $idx = "All" }

        #I really want to get right to the transcript and peptide pages, but
        #can't see how to do that without knowing the org name too, which
        #I don't know at this point. (ensembl org name, not real org name)

        return "<a href='http://www.ensembl.org/Homo_sapiens/searchview?species=all&idx=$idx&q=$key' target=_blank>$ensembl</a>";
    }
    return $ensembl;
}

sub entrezgene_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$entrezgene) = @_;

    if ($entrezgene =~ /^EntrezGene:(\S+)$/)
    {
        return "<a href='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene&cmd=Retrieve&dopt=full_report&list_uids=$1' target=_blank>$entrezgene</a>";
    }
    return $entrezgene;
}

sub mim_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$mim) = @_;

    if ($mim =~ /^MIM:(\S+)$/)
    {
        return "<a href='http://www3.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=$1' target=_blank>$mim</a>";
    }
    return $mim;
}

sub hgnc_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$hgnc) = @_;

    if ($hgnc =~ /^HGNC:(\S+)$/)
    {
        return "<a href='http://www.gene.ucl.ac.uk/cgi-bin/nomenclature/searchgenes.pl?field=symbol&anchor=equals&match=$1&symbol_search=Search&number=50&format=html&sortby=symbol' target=_blank>$hgnc</a>";
    }

    return $hgnc;
}

sub unigene_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$unigene) = @_;

    if ($unigene =~ /^UniGene:(\S+)$/)
    {
        return "<a href='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=unigene&cmd=search&term=$1' target=_blank>$unigene</a>";
    }
    return $unigene;
}

sub ipi_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$ipi) = @_;

    if ($ipi =~ /^IPI:(\S+)$/)
    {
        return "<a href='http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-id+AEoS1R8Jnn+-e+[IPI:\'$1\']+-qnum+1+-enum+1' target=_blank>$ipi</a>";
    }
    return $ipi;
}

sub wp_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$wp) = @_;

    #wormbase

    if ($wp =~ /^WP:(\S+)$/)
    {
        return "<a href='http://www.wormbase.org/db/searches/basic?class=Any&query=$1&Search=Search' target=_blank>$wp</a>";
    }
    return $wp;
}

sub fb_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$fb) = @_;

    #flybase

    if ($fb =~ /^FB:(\S+)$/)
    {
        return "<a href='http://flybase.bio.indiana.edu/.bin/fbidq.html?$1' target=_blank>$fb</a>";
    }
    return $fb;
}

sub fborf_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$fb) = @_;

    #flybase

    if ($fb =~ /^FlyBaseORFNames:(\S+)$/)
    {
        return "<a href='http://flybase.bio.indiana.edu/.bin/fbidq.html?$1' target=_blank>$fb</a>";
    }
    return $fb;
}

sub sgd_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$sgd) = @_;

    #yeast

    if ($sgd =~ /^SGD_LOCUS:(\S+)$/)
    {
        return "<a href='http://db.yeastgenome.org/cgi-bin/locus.pl?locus=$1' target=_blank>$sgd</a>";
    }
    return $sgd;
}




sub set_map_links {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$x) = @_;
    my($before,$match,$after);

    my $org = ($cgi->param('org') || $cgi->param('genome') || "");

    if ($x =~ /^(.*)(MAP\d+)(.*)/s)
    {
        $before = $1;
        $match = $2;
        $after = $3;
        return &set_map_links($cgi,$before) . &map_link($cgi,$match,$org) . &set_map_links($cgi,$after);
    }
    return $x;
}



sub map_link {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cgi,$map,$org) = @_;

    my $user = $cgi->param('user');
    $user = $user ? $user : "";
    $org = $org ? $org : "";

    my $url = "show_kegg_map.cgi?user=$user&map=$map&org=$org";
#rel    my $url = &FIG::cgi_url() . "/show_kegg_map.cgi?user=$user&map=$map&org=$org";
    my $link = "<a href=\"$url\">$map</a>";
    return $link;
}


sub java_buttons {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    ## ADDED BY RAE
    # Provides code to include check all/first half/second half/none for javascrspt
    # this takes two variables - the form name provided in start_form with the
    # -name => field and the checkbox name
    my ($form, $button)=@_;

    return  qq(<input type="button" name="CheckAll" value="Check All"\nonClick="checkAll(document.$form.$button)">\n)
          . qq(<input type="button" name="CheckFirst" value="Check First Half"\nonClick="checkFirst(document.$form.$button)">\n)
          . qq(<input type="button" name="CheckSecond" value="Check Second Half"\nonClick="checkSecond(document.$form.$button)">\n)
          . qq(<input type="button" name="UnCheckAll" value="Uncheck All"\nonClick="uncheckAll(document.$form.$button)">\n);
}


#
#  Provides code to include check all/first half/second half/none for javascrspt
#  this takes two variables - the form name provided in start_form with the
#  -name => field and the checkbox name.
#
sub java_buttons_ext {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my ( $form, $button ) = @_;

    return  qq(<input type="button" name="CheckAll"    value="Check All"             onClick="checkAll(document.$form.$button)">\n)
          . qq(<input type="button" name="CheckFirst"  value="Check First Half"      onClick="checkFirst(document.$form.$button)">\n)
          . qq(<input type="button" name="CheckSecond" value="Check Second Half"     onClick="checkSecond(document.$form.$button)">\n)
          . qq(<input type="button" name="CheckToMark" value="Check to First Marked" onClick="checkToFirstMarked1(document.$form.$button)">\n)
          . qq(<input type="button" name="UnCheckAll"  value="Uncheck All"           onClick="uncheckAll(document.$form.$button)">\n);
}

=head3 sub_link

    my $htmlText = HTML::sub_link($cgi, $sub, gid);

Create a subsystem link. The link will be to the display page if there is no
user or we are in SPROUT mode; otherwise it will be to the edit page.

=over 4

=item cgi

CGI query object for the current web session. The parameters of special interest
are C<SPROUT> and C<user>. If the user is non-blank and SPROUT mode is 0, then
the subsystem's edit page will be shown rather than its display page.

=item sub

Name of the desired subsystem. It will be cleaned of underscores before the
hyperlink is applied.

=item gid

Genome ID to be specified as the focus.

=back

=cut

sub sub_link {
    # Allow call as an instance in addition to the authorized method.
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    # Get the parameters.
    my ($cgi, $sub, $gid) = @_;
    # Declare the return variable.
    my $retVal;
    # Clean the subsystem name for display purposes. This is a very
    # different thing from URL-escaping.
    my $cleaned = Tracer::Clean($sub);
    $cleaned =~ s/_/ /g;
    # URL-escape the subsystem name for use in the link.
    my $linkable = uri_escape($sub);
    # Determine the mode. Note we use the little OR trick to insure that
    # we have the correct value for plugging into the output link.
    my $user = $cgi->param('user') || "";
    my $sproutMode = $cgi->param('SPROUT') || 0;
    if ($user && ! $sproutMode) {
        # A SEED user is calling, so we go to the edit page.
        $retVal = "<a href=\"subsys.cgi?ssa_name=$linkable&request=show_ssa&user=$user&can_alter=1&check=1&sort=&show_clusters=&show_minus1=1\">$cleaned</a>";
    } else {
        # A visitor or SPROUT user is calling, so we go to the display page.
        $retVal = "<a href=\"$FIG_Config::cgi_url/seedviewer.cgi?page=Subsystems;subsystem=$linkable\">$cleaned</a>";
    }
    # Return the result.
    return $retVal;
}

sub reaction_map_link {
    my($mapID, @reaction_list) = @_;
    if($mapID =~ /\d+/)
    {
        my $reactions = join "+", @reaction_list;
        if ($reactions ne "")
        {
            $reactions = "+".$reactions;
        }

        return "<a href=http://www.genome.jp/dbget-bin/show_pathway?rn$mapID$reactions>$mapID</a>";
    }
    else
    {
        return $mapID;
    }
}

sub compound_link {
    my($compound) = @_;
    if($compound =~ /^C\d+/)
    {
        return "<a href=\"javascript:void(0)\"onclick=\"window.open('http://www.genome.jp/dbget-bin/www_bget?compound+$compound','$&','height=640,width=800,scrollbars=yes,toolbar=yes,status=yes,resizable=yes')\">$compound</a>";
    }
    else
    {
        return $compound;
    }
}


sub reaction_link {
    my($reaction) = @_;
    if ($reaction =~ /^(\*)?(R\d+)/)
    {
        # return "<a href=\"http://www.genome.ad.jp/dbget-bin/www_bget?rn+$2\" target=reaction$$>$reaction</a>";
        return "<a href=\"javascript:void(0)\"onclick=\"window.open('http://www.genome.ad.jp/dbget-bin/www_bget?rn+$reaction','$&','height=640,width=800,scrollbars=yes,toolbar=yes,status=yes,resizable=yes')\">$reaction</a>";
    }
    return $reaction;
}


sub html_for_assignments {
    my($fig,$user,$peg_sets) = @_;
    my $i;

    my @vals = ();
    my $set = 1;
    foreach my $peg_set (@$peg_sets)
    {
        for ($i=0; ($i < @$peg_set); $i++)
        {
            my $peg = $peg_set->[$i];
            push(@vals,'show=' . join("@",($set,$i+1,$peg,&FIG::abbrev($fig->org_of($peg)),"")));
        }
        $set++;
    }

    $ENV{'REQUEST_METHOD'} = 'GET';
    $ENV{'QUERY_STRING'} = join('&',('request=show_commentary',"uni=1","user=$user",@vals));
    my $out = join("",`$FIG_Config::fig/CGI/chromosomal_clusters.cgi`);
    $out =~ s/^.*?<form/<form/si;
    $out =~ s/^(.*)<table.*/$1/si;
    return $out;
}

=head1 rss_feed

Add something to the RSS feed. The rss feeds are stored in the Html directory, and there are several RSS feeds:
        SEED.rss                - everything gets written here
        SEEDgenomes.rss                 - whenever a genome is added to the SEED
        SEEDsubsystems.rss      - whenever a subsystem is edited (or should this be added?)


RSS feeds must contain a title, description, and link. The title is what is seen e.g. from the firefox or safari pull down menu. The description is seen from within an rss aggregator, and may be displayed on web pages and so on.

The method takes a reference to an array containing the file names for the RSS feeds to add your item to, and a hash of items for the xml. Only title, description, and link are required tags in the XML.

The file names are the full name of the file, eg SEEDsubsystems.rss, SEEDgenomes.rss. Be aware that this is a file name, though, so don't uses special characters. The path will be added.

The has can have these keys:

REQUIRED:
title       : the title. This is usually what is seen by the user in the pull down menu
description : a more complete description that is often seen is rss viewers but not always
link        : link to the item that was added/edited
All other keys are treated as optional RSS arguments and written to the file.

At most, $max_entries recent entries are stored in the rss file, and this is currently 50.

RSS files are quite simple, and contain some standard header information, and then individual items surrounded by an <item> </item> tag. Note that there is also an initial title/description/link set that describes the file.


=cut

sub rss_feed {
 shift if UNIVERSAL::isa($_[0],__PACKAGE__);
 my ($files, $args)=@_;

 # how many entries to store in the file
 my $max_entries=50;

 foreach my $a (keys %$args) {if ($a =~ /^-(.*)/) {my $b=$1; $args->{$b}=$args->{$a}; delete $args->{$a}}}

 my $filepath=$FIG_Config::fig."/CGI/Html/rss";
 # check for the directory and if not, make it
 mkdir $filepath unless (-d $filepath);

 # note that $info is a hash of references to hashes that are written out as headers in the file
 my $info=
 {
  "SEED.rss" =>
   {
        title           => "The SEED",
        description     => "Latest news from the SEED",
        link            => "Html/rss/SEED.rss",
   },

  "SEEDsubsystems.rss" =>
  {
        title           => "SEED Subsystems",
        description     => "Recently updated SEED subsystems",
        link            => "Html/rss/SEEDsubsystems.rss",
  },

  "SEEDsubsystems.rss" =>
  {
        title           => "SEED Genomes",
        description     => "Genomes recently added to the SEED",
        link            => &FIG::cgi_url()."/Html/rss/SEEDsubsystems.rss",
  },

 };


 # build the new xml
 my $xml = "\t<item>\n";
 foreach my $qw ("title", "description", "link") {
  unless ($args->{$qw}) {
   print STDERR "You need to include a $qw tag in your RSS description\n";
   return(0);
  }
  # we need to do something a bit funky with the link. We can't have ampersands in the <link> </link> in valid html
  # so we are going to pull out the links and uri_escape just the part after the .cgi
  if ($qw eq "link")
  {
   $args->{$qw} =~ /^(.*?\.cgi.)(.*)$/;
   $args->{$qw} = $1.uri_escape($2) if ($1 && $2);
  }

  $xml .= "\t\t<$qw>".$args->{$qw}."</$qw>\n";
  delete $args->{$qw};
 }

 foreach my $tag (grep {!/type/i} keys %$args)
 {
  $xml .= "\t\t<$tag>".$args->{$tag}."</$tag>\n";
 }

 $xml .= "\t</item>\n";


 my @files=("SEED.rss");
 if ($args->{"type"}) {
    my $type = $args->{type};
    push @files, "SEED.$type.rss"
}

 foreach my $file ("SEED.rss", @$files)
 {
  if (-e "$filepath/$file")
  {
   my @out; # the new content of the file
   my $itemcount=0; # how many <item> </item>'s are we keeping
   my $initem; # are we in an item?
   open(IN, "$filepath/$file") || die "Can't open $filepath/$file";
   while (<IN>)
   {
    if (/\<item\>/) {
     push @out, $xml, unless ($itemcount);
     $itemcount++;
     $initem=1;
    }
    if (/\<\/item\>/) {$initem=0; next if ($itemcount > $max_entries)}
    next if ($initem && $itemcount > $max_entries);
    push @out, $_;
   }
   close IN;
   open(OUT, ">$filepath/$file") || die "Can't open $filepath/$file for writing";
   print OUT @out;
  }
  else
  {
   open(OUT, ">$filepath/$file") || die "Can't open $filepath/$file for writing";
   print OUT "<?xml version=\"1.0\"?>\n<rss version=\"2.0\">\n<channel>\n";
   if ($info->{$file})
   {
     # we're going to sanity check each of the three options we output, just to be sure
     foreach my $qw ("title", "description", "link")
     {
       if ($info->{$file}->{$qw})
       {
          print OUT "<$qw>", $info->{$file}->{$qw}, "</$qw>\n";
       } else {
          print STDERR "Please add a $qw for $file\n"; print OUT "<$qw>$file</$qw>\n";
       }
     }
   }
   else {
    print STDERR "Please define title, link, and description information for $file\n";
    print OUT "<title>$file</title>\n<description>An RSS feed</description>\n<link>", &FIG::cgi_url, "</link>\n";
   }
   print OUT "\n", $xml;
   print OUT "\n", "</channel>\n</rss>\n"
  }
 }
}



1;

