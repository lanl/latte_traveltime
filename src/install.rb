#!/usr/bin/env ruby

for dir in ["2d", "3d"]

        if (not ARGV.empty?) and ARGV[0] == "clean"
                system(" cd #{dir}; make clean; make &")
        else
                system(" cd #{dir}; make & ")
        end

end
