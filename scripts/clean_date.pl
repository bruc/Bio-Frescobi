#!perl

use strict;
use warnings;

# Filter to mask dates:

while (<STDIN>) {
    s/(Mon|Tue|Wed|Thu|Fri|Sat|Sun) (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [ 0-9]{1,2} \d{2}:\d{2}:\d{2} [A-Z]{3,5} \d{4}/xxx xxx xx xx:xx:xx XXX xxxx/g;
    s/(Mon|Tue|Wed|Thu|Fri|Sat|Sun) (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [ 0-9]{1,2} \d{2}:\d{2}:\d{2} \d{4} [A-Z]{3,5}/xxx xxx xx xx:xx:xx XXX xxxx/g;
    s/\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}(.\d+)?-\d{2}/xxxx-xx-xx xx:xx:xx-xx/g;
    print $_;
}
