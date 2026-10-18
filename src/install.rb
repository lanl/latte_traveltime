#!/usr/bin/env ruby

# Usage: ruby install.rb [clean] [name=value ...]
#
# clean removes the old objects and executables before compiling. Each
# name=value goes to make and overrides the variable of the same name in the
# Makefiles, e.g.,
#
#   ruby install.rb clean use_cuda=on cuda_arch=sm_80
#
# also compiles the CUDA executables of the 3D solvers for the GPU
# architecture sm_80.

require 'shellwords'

clean = ARGV.include?("clean")
options = ARGV.grep(/\A\w+=/)
unknown = ARGV - ["clean"] - options
if not unknown.empty?
        abort("Unknown argument: #{unknown.join(' ')}\nUsage: ruby install.rb [clean] [name=value ...]")
end
opts = Shellwords.join(options)

for dir in ["2d", "3d"]

        if clean
                system(" cd #{dir}; make clean #{opts}; make #{opts} &")
        else
                system(" cd #{dir}; make #{opts} & ")
        end

end
