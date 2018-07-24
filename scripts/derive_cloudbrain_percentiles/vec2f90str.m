function [ string ] = vec2f90str(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
  n = length (x);
  string = sprintf ('%e',x(1));
  for k=2:n
      string = sprintf ('%s, %e ',string,x(k));
  end

end
