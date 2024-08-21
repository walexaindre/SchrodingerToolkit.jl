#What a report should contain ?
#1. Description of method and parameters to be tested
#2. Some plots of the results
#3. A table of the results
#4. A conclusion


# How to do this ?
# I need to define a list of mesh sizes and a fixed time step size


function scientific_notation_to_latex(str::String;inline_math::Bool = true)
    comp = split(str, "e")
    left = comp[1]
    right = comp[2]
    idx = 1
    #Deleting leading zeros
    for elem in right[2:end-1]
        if elem == '0'
            idx += 1
            continue
        else
            break
        end
    end

    rightnum = right[idx+1:end]
    sign = right[1]

    if sign == '+'
        sign = ""
    end

    innerlatex = "$left \\times 10^{$(sign)$rightnum}"

    if inline_math
        return "\$$innerlatex\$"
    else
        return innerlatex
    end
    error("Unreachable code")
end


function scientific_notation_to_latex(int::Number,inline_math::Bool = true)
    str = @sprintf("%.4e",int)
    return scientific_notation_to_latex(str,inline_math = inline_math) 
end

function Table(header::Vector{String},data::Array{String,2})
    columns = length(header)
    rows = size(data,1)
    if columns != size(data,2)
        error("Number of columns in header and data do not match")
    end
    table_prelude = "\\begin{table}[H]\n\\centering\n\\begin{tabular}{"
    table_prelude *= "l"
    for i in 2:columns
        table_prelude *= "c"
    end
    table_prelude *= "}\n"
    table_header = "\\hline"*join(table_header, " & ")*"\\\\\n\\hline\\hline\n"

    table_data = ""

    for i in 1:rows
        table_data*=join(data[i,:], " & ")*"\\\\\n"
    end

    table_end = "\\hline\n\\end{tabular}\n\\end{table}"

    return table_prelude*table_header*table_data*table_end
end
