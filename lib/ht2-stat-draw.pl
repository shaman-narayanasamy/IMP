#!/usr/bin/perl

#     HTQC - a high-throughput sequencing quality control toolkit
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use File::Path;
use File::Spec::Functions;
use Getopt::Long;

# output file name base
my $CYC_QUAL_1       = 'cycle_quality_1';
my $CYC_QUAL_2       = 'cycle_quality_2';
my $CYC_QUAL_1_BOX   = 'cycle_quality_box_1';
my $CYC_QUAL_2_BOX   = 'cycle_quality_box_2';
my $READS_QUAL       = 'reads_quality';
my $LANE_TILE_QUAL_1 = 'lane_tile_quality_1';
my $LANE_TILE_QUAL_2 = 'lane_tile_quality_2';
my $CYC_COMP_1       = 'cycle_composition_1';
my $CYC_COMP_2       = 'cycle_composition_2';
my $UMASK_LEN        = 'reads_length';
my $QUAL_QQ          = 'quality_QQ';

my $dir_data;
my $WIDTH    = 600;
my $HEIGHT   = 400;
my $SZ_TITLE = 12;
my $SZ_TEXT  = 9;
my $FORMAT   = 'png';

GetOptions(
    'dir=s'        => \$dir_data,
    'format=s'     => \$FORMAT,
    'width=s'      => \$WIDTH,
    'height=s'     => \$HEIGHT,
    'title-size=s' => \$SZ_TITLE,
    'text-size=s'  => \$SZ_TEXT,
    'help'         => \&show_help
);

#
# validate dir
#
die "data directory not specified" if !defined $dir_data;
die "data directory not exist" if !-d $dir_data;

my $file_info = catfile $dir_data, 'info.tab';
die "info file not exist" if !-f $file_info;

#
# get info
#
open INFO, '<', $file_info or die "failed to open info file '$file_info': $!";
my $PAIRED;
my $LENGTH;
my $QUAL_FROM;
my $QUAL_TO;
my $MASK;

while (<INFO>) {
    chomp;
    my ( $key, $value ) = split /\t/;
    if ( $key eq 'paired' ) {
        $PAIRED = $value eq 'T' ? 1 : 0;
    }
    elsif ( $key eq 'length' ) {
        $LENGTH = $value;
    }
    elsif ( $key eq 'quality range' ) {
        $value =~ /^(\d+)-(\d+)$/
          or die "failed to parse quality range: '$value'";
        $QUAL_FROM = $1;
        $QUAL_TO   = $2;
    }
    elsif ( $key eq 'mask' ) {
        $MASK = $value eq 'T' ? 1 : 0;
    }
}

close INFO;

#
# test the existence of gnuplot
#
my $HAS_GNUPLOT = 1;
`gnuplot -V`;
if ( $? != 0 ) {
    print STDERR <<INFO_GNUPLOT_NOTFOUND;
Gnuplot is not found in your system. Plot scripts are generated but are not
rendered.
Please install Gnuplot using the package management system of your Unix
distribution, or visit Gnuplot website: http://sourceforge.net/projects/gnuplot
INFO_GNUPLOT_NOTFOUND
    $HAS_GNUPLOT = 0;
}

#
# write gnuplot scripts
#
plot_cycle_quality_heatmap( $CYC_QUAL_1, 1 );
plot_cycle_quality_heatmap( $CYC_QUAL_2, 2 ) if $PAIRED;

plot_cycle_quality_box( $CYC_QUAL_1_BOX, 1 );
plot_cycle_quality_box( $CYC_QUAL_2_BOX, 2 ) if $PAIRED;

plot_reads_quality($READS_QUAL);

plot_lane_tile_quality( $LANE_TILE_QUAL_1, 1 );
plot_lane_tile_quality( $LANE_TILE_QUAL_2, 2 ) if $PAIRED;

plot_cycle_composition( $CYC_COMP_1, 1 );
plot_cycle_composition( $CYC_COMP_2, 2 ) if $PAIRED;

plot_reads_length($UMASK_LEN);

plot_quality_qq($QUAL_QQ) if $PAIRED;

#
# subs
#

# gnuplot scripts
sub plot_cycle_quality_heatmap {
    my $base = shift;
    my $side = shift;

    my $file_tab   = catfile $dir_data, $base . '.tab';
    my $file_plot  = catfile $dir_data, $base . '.gnuplot';
    my $file_chart = catfile $dir_data, $base . '.' . $FORMAT;
    my $file_dat   = catfile $dir_data, $base . '.dat';

    # read table
    # write data file for gnuplot
    open my $h_tab, '<', $file_tab
      or die "failed to open '$file_tab' for read: $!";
    open my $h_dat, '>', $file_dat
      or die "failed to open '$file_dat' for write: $!";

    my @qual_list;

    my $line = 0;
    while (<$h_tab>) {
        chomp;
        if ( $line == 0 ) {
            ( undef, @qual_list ) = split /\t/;
        }
        else {
            my ( $cyc, @values ) = split /\t/;
            for ( my $i = 0 ; $i < @values ; $i++ ) {
                print $h_dat join( "\t", $cyc, $qual_list[$i], $values[$i] ),
                  "\n";
            }
        }
    }
    continue { $line++ }

    close $h_tab;
    close $h_dat;

    # write script
    my $x_min = -0.5;
    my $x_max = $LENGTH + 0.5;
    my $y_min = $QUAL_FROM - 0.5;
    my $y_max = $QUAL_TO + 0.5;

    open my $h_plot, '>', $file_plot or die $!;
    print $h_plot <<HEREDOC;
set title "Cycle Quality of Reads $side" font "*,$SZ_TITLE"
unset key
set terminal $FORMAT size $WIDTH,$HEIGHT
set output "$file_chart"
set palette rgbformulae 21,22,23

set cblabel "Num reads" font "*,$SZ_TEXT"
set cbtics font "*,$SZ_TEXT"

set xrange [$x_min:$x_max]
set yrange [$y_min:$y_max]
set xlabel "Cycle" font "*,$SZ_TEXT"
set ylabel "Quality" font "*,$SZ_TEXT"
set xtics font "*,$SZ_TEXT"
set ytics font "*,$SZ_TEXT"

plot '$file_dat' using 1:2:3 with image

HEREDOC
    close $h_plot;

    # run gnuplot
    if ($HAS_GNUPLOT) {
        system( 'gnuplot', $file_plot ) == 0 or die "gnuplot failed";
    }
}

sub plot_cycle_quality_box {
    my $base = shift;
    my $side = shift;

    my $file_tab   = catfile $dir_data, $base . '.tab';
    my $file_plot  = catfile $dir_data, $base . '.gnuplot';
    my $file_chart = catfile $dir_data, $base . '.' . $FORMAT;

    # write script
    my $x_min = -0.5;
    my $x_max = $LENGTH + 0.5;
    my $y_min = $QUAL_FROM - 0.5;
    my $y_max = $QUAL_TO + 0.5;

    open my $h_plot, '>', $file_plot or die $!;
    print $h_plot <<HEREDOC;
set title "Cycle Quality of Reads $side" font "*,$SZ_TITLE"
unset key
set terminal $FORMAT size $WIDTH,$HEIGHT
set output "$file_chart"

set xrange [$x_min:$x_max]
set yrange [$y_min:$y_max]
set xlabel "Cycle" font "*,$SZ_TEXT"
set ylabel "Quality" font "*,$SZ_TEXT"
set xtics font "*,$SZ_TEXT"
set ytics font "*,$SZ_TEXT"

set boxwidth 0.8 relative

plot '$file_tab' using 1:3:2:6:5 with candlesticks lc rgbcolor '#00007f',\\
     '$file_tab' using 1:7 with lines lc rgbcolor 'red'

HEREDOC
    close $h_plot;

    # run gnuplot
    if ($HAS_GNUPLOT) {
        system( 'gnuplot', $file_plot ) == 0 or die "gnuplot failed";
    }
}

sub plot_reads_quality {
    my $base = shift;

    my $file_tab   = catfile $dir_data, $base . '.tab';
    my $file_plot  = catfile $dir_data, $base . '.gnuplot';
    my $file_chart = catfile $dir_data, $base . '.' . $FORMAT;

    # write script
    open my $h_plot, '>', $file_plot or die $!;
    print $h_plot <<HEREDOC;
set title "Reads Quality" font "*,$SZ_TITLE"
set terminal $FORMAT size $WIDTH,$HEIGHT

set key font "*,$SZ_TEXT"

set xlabel "Quality" font "*,$SZ_TEXT"
set ylabel "Num reads" font "*,$SZ_TEXT"
set xtics font "*,$SZ_TEXT"
set ytics font "*,$SZ_TEXT"

set output "$file_chart"

HEREDOC

    print $h_plot
"plot '$file_tab' using 1:2 lc rgbcolor '#ff0000' title 'read 1 quality' with lines,\\\n";
    print $h_plot
"     '$file_tab' using 1:3 lc rgbcolor '#7f0000' title 'read 1 quality accum' with lines";

    if ($PAIRED) {
        print $h_plot ",\\\n";
        print $h_plot
"     '$file_tab' using 1:4 lc rgbcolor '#0000ff' title 'read 2 quality' with lines,\\\n";
        print $h_plot
"     '$file_tab' using 1:5 lc rgbcolor '#00007f' title 'read 2 quality accum' with lines\n";
    }
    else {
        print $h_plot "\n";
    }

    close $h_plot;

    # run gnuplot
    if ($HAS_GNUPLOT) {
        system( 'gnuplot', $file_plot ) == 0 or die "gnuplot failed";
    }
}

sub plot_lane_tile_quality {
    my $base = shift;
    my $side = shift;

    my $file_tab   = catfile $dir_data, $base . '.tab';
    my $file_plot  = catfile $dir_data, $base . '.gnuplot';
    my $file_chart = catfile $dir_data, $base . '.' . $FORMAT;
    my $file_dat   = catfile $dir_data, $base . '.dat';

    # read ticks
    # create dat file
    my @titles;
    my @ticks_list;

    open my $h_tab, '<', $file_tab
      or die "failed to open '$file_tab' for read: $!";
    open my $h_dat, '>', $file_dat
      or die "failed to open '$file_dat' for write: $!";
    my $n = 0;
    while (<$h_tab>) {
        chomp;
        if ( $n == 0 ) {
            ( undef, @titles ) = split /\t/;
        }
        else {
            my ( $lane, $tile, @data ) = split /\t/;
            push @ticks_list, "\"lane $lane tile $tile\" $n";
            print $h_dat join( "\t", $n, @data ), "\n";
        }
    }
    continue { $n++ }
    close $h_tab;
    close $h_dat;

    my $ticks_str = join ', ', @ticks_list;

    # write script
    open my $h_plot, '>', $file_plot or die $!;
    print $h_plot <<HEREDOC;
set title "Lane-Tile Quality for Reads $side" font "*,$SZ_TITLE"
set terminal $FORMAT size $WIDTH,$HEIGHT

set key outside font "*,$SZ_TEXT"

set xlabel "Tiles" font "*,$SZ_TEXT"
set ylabel "Num reads" font "*,$SZ_TEXT"
set xtics rotate 90 ($ticks_str) font "*,$SZ_TEXT"
set ytics font "*,$SZ_TEXT"

set style fill solid
set output "$file_chart"
plot '$file_dat' using 1:(\$5+\$4+\$3+\$2):(1) with boxes lc rgbcolor "blue" title "quality > 30", \\
     '$file_dat' using 1:(\$4+\$3+\$2):(1) with boxes lc rgbcolor "green" title "20 < quality < 30", \\
     '$file_dat' using 1:(\$3+\$2):(1) with boxes lc rgbcolor "red" title "10 < quality < 20", \\
     '$file_dat' using 1:2:(1) with boxes lc rgbcolor "black" title "quality < 10"
HEREDOC

    close $h_plot;

    # run gnuplot
    if ($HAS_GNUPLOT) {
        system( 'gnuplot', $file_plot ) == 0 or die "gnuplot failed";
    }
}

sub plot_cycle_composition {
    my $base = shift;
    my $side = shift;

    my $file_tab   = catfile $dir_data, $base . '.tab';
    my $file_plot  = catfile $dir_data, $base . '.gnuplot';
    my $file_chart = catfile $dir_data, $base . '.' . $FORMAT;

    my $x_min = -0.5;
    my $x_max = $LENGTH + 0.5;

    open my $h_plot, '>', $file_plot or die $!;
    print $h_plot <<heredoc;
set title "Cycle Composition for Read $side" font "*,$SZ_TITLE"
set terminal $FORMAT size $WIDTH,$HEIGHT
set style fill solid

set key outside font "*,$SZ_TEXT"

set xrange [$x_min:$x_max]
set xlabel "Cycle" font "*,$SZ_TEXT"
set ylabel "Num reads" font "*,$SZ_TEXT"
set xtics font "*,$SZ_TEXT"
set ytics font "*,$SZ_TEXT"

set output "$file_chart"
heredoc
    print $h_plot "plot";

    if ($MASK) {
        print $h_plot <<heredoc;
'$file_tab' using 1:(\$7+\$6+\$5+\$4+\$3+\$2) with boxes lc rgbcolor "yellow" title "MASK", \\
heredoc
    }

    print $h_plot <<heredoc;
'$file_tab' using 1:(\$6+\$5+\$4+\$3+\$2) with boxes lc rgbcolor "grey" title "N", \\
     '$file_tab' using 1:(\$5+\$4+\$3+\$2) with boxes lc rgbcolor "black" title "C", \\
     '$file_tab' using 1:(\$4+\$3+\$2) with boxes lc rgbcolor "blue" title "G", \\
     '$file_tab' using 1:(\$3+\$2) with boxes lc rgbcolor "red" title "T", \\
     '$file_tab' using 1:2 with boxes lc rgbcolor "green" title "A"
heredoc
    close $h_plot;

    # run gnuplot
    if ($HAS_GNUPLOT) {
        system( 'gnuplot', $file_plot ) == 0 or die "gnuplot failed";
    }
}

sub plot_reads_length {
    my $base = shift;

    my $file_tab   = catfile $dir_data, $base . '.tab';
    my $file_plot  = catfile $dir_data, $base . '.gnuplot';
    my $file_chart = catfile $dir_data, $base . '.' . $FORMAT;

    open my $h_plot, '>', $file_plot or die $!;
    print $h_plot <<heredoc;
set title "Reads Length" font "*,$SZ_TITLE"
set terminal $FORMAT size $WIDTH,$HEIGHT

set key font "*,$SZ_TEXT"

set xlabel "Length" font "*,$SZ_TEXT"
set ylabel "Num reads" font "*,$SZ_TEXT"
set xtics font "*,$SZ_TEXT"
set ytics font "*,$SZ_TEXT"

set output "$file_chart"
heredoc
    print $h_plot
"plot '$file_tab' using 1:2  lc rgbcolor '#ff0000' title 'read 1 length' with lines,\\\n";
    print $h_plot
"     '$file_tab' using 1:3  lc rgbcolor '#7f0000' title 'read 1 length accum' with lines";
    if ($PAIRED) {
        print $h_plot ",\\\n";
        print $h_plot
"     '$file_tab' using 1:4 lc rgbcolor '#0000ff' title 'read 2 length' with lines,\\\n";
        print $h_plot
"     '$file_tab' using 1:5 lc rgbcolor '#00007f' title 'read 2 length accum' with lines\n";
    }
    else {
        print $h_plot "\n";
    }
    close $h_plot;

    # run gnuplot
    if ($HAS_GNUPLOT) {
        system( 'gnuplot', $file_plot ) == 0 or die "gnuplot failed";
    }
}

sub plot_quality_qq {
    my $base = shift;

    my $file_tab   = catfile $dir_data, $base . '.tab';
    my $file_plot  = catfile $dir_data, $base . '.gnuplot';
    my $file_chart = catfile $dir_data, $base . '.' . $FORMAT;

    my $label_x = $QUAL_FROM + 3;
    my $label_y = $QUAL_TO - 3;

    # get pearson r value
    open my $h_tab, '<', $file_tab
      or die "failed to open '$file_tab' for read: $!";
    my $tab_header = <$h_tab>;
    chomp $tab_header;
    $tab_header =~ /# pearson correlation: (\S+)/
      or die "failed to parse pearson correlation from '$tab_header'";
    my $pearson_text = $1;
    close $h_tab;

    # write gnuplot script
    open my $h_plot, '>', $file_plot or die $!;
    print $h_plot <<heredoc;
set title "Read Pair Quality QQ plot" font "*,$SZ_TITLE"
set terminal $FORMAT size $WIDTH,$HEIGHT
set output "$file_chart"
unset key

set xrange [$QUAL_FROM:$QUAL_TO]
set yrange [$QUAL_FROM:$QUAL_TO]
set xtics font "*,$SZ_TEXT"
set ytics font "*,$SZ_TEXT"

set label 'r = $pearson_text' at $label_x,$label_y font "*,$SZ_TEXT"
plot '$file_tab' with lines

heredoc
    close $h_plot;

    # run gnuplot
    if ($HAS_GNUPLOT) {
        system( 'gnuplot', $file_plot ) == 0 or die "gnuplot failed";
    }
}

# help document
sub show_help {
    print <<heredoc;
$0 - render ht_stat outputs

Usage:
  ht_stat_draw.pl --dir STAT_DIR

Options:
--dir         the output directory of ht_stat.
--format      picture format. [$FORMAT]
--width       picture width. [$WIDTH]
--height      picture height. [$HEIGHT]
--title-size  font size for picture titles. [$SZ_TITLE]
--text-size   font size for axis ticks, legends, etc.. [$SZ_TEXT]
--help        show help
heredoc
    exit(0);
}
