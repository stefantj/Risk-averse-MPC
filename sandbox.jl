using Convex


##################
# form_variables #
##################
# recursively declares variables
# depth   - current recursion depth
# m       - total number of levels
# varname - base variable name to instantiate
# vsize1  - 1st dimension of variable
# vsize2  - 2nd dimension of variable
function form_variables(depth, m, varname, vsize1, vsize2)
    depth -= 1;
    for i = 1:m
       #Extend variable subindex
       vname = string(varname, "_",i)

       # Form expression:
       v = parse(string(vname,"=Variable(",vsize1,",",vsize2,")"))
       eval(v);
       if(depth > 0)
          form_variables(depth, m, vname, vsize1, vsize2);
       end
    end
end


#######################
# BEGIN SCRIPT        #
#######################

depth = 3;
m = 10;
numvars = 0;
for k = 1:depth
   numvars += m^k;
end

#Forms variables Y_1:m, Y_1:m_1:m, Y_1:m_1:m_1:m
println("Forming ", numvars, " variables")
@time form_variables(depth, m, "global Y", 10, 10);

W_N = Variable(m,m);

print("Forming constraints")


