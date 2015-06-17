using Convex;


# form_variables(depth, m, varname, vsize1, vsize2)
#
# Recursively declares variables
# depth   - current recursion depth
# m       - Number of variables per level
# varname - base variable name to increment on
# vsize1  - 1st dimension of variable
# vsize2  - 2nd dimension of variable
function form_variables(depth, m, varname, vsize1, vsize2)
   depth-=1;
   for i = 1:m
      vname = string(varname, "_", i);
      v = parse(string(vname, "=Variable(", vsize1, ",", vsize2, ");"));
      eval(v);
      if(depth > 0)
          form_variables(depth, m, vname, vsize1, vsize2);
      end
   end
end


