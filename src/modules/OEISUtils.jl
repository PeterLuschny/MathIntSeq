# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module OEISUtils
using Requests, URIParser, SeqBase, Nemo

export oeis_writebfile, oeis_trimdata, oeis_remote, oeis_local, oeis_isinstalled
export oeis_notinstalled, oeis_path, oeis_search

#doc"""
#Directory of oeis data.
#"""
MODDIR = realpath(joinpath(dirname(@__FILE__)))
DATADIR = joinpath(dirname(MODDIR), "data")

doc"""
Returns the path where the oeis data is expected.
"""
oeis_path() = joinpath(DATADIR, "stripped")

doc"""
Indicates if the local copy of the OEIS data (the so-called
'stripped' file) is installed (in MathIntSeq/data).
"""
oeis_isinstalled() = isfile(oeis_path())

doc"""
Indicates if the local copy of the OEIS data (the so-called
'stripped' file) is not installed and warns.
"""
function oeis_notinstalled()
    if !oeis_isinstalled()
        warn("OEIS data not installed! Download stripped.gz from oeis.org,")
        warn("expand it and put it in the directory MathIntSeq/data.")
        return true
    end
    return false
end

doc"""
Write a so-called b-file for submission to the OEIS. The file is saved in the
'data' directory.
"""
function oeis_writebfile(anum, fun, offset::Int, len::Int)

    if !ismatch(r"^A[0-9]{6}$", anum)
        warn("Not a valid A-number!")
        return
    end

    filename = joinpath(DATADIR, "b" * anum[2:end] * ".txt")
    info("Writing " * anum * " to " * filename)

    f = open(filename, "w")
    for n in offset:(offset + len)
        println(f, n, " ", fun(n))
    end
    close(f)
end

function oeis_writebfile(anum, list)

    if !ismatch(r"^A[0-9]{6}$", anum)
        warn("Not a valid A-number!")
        return
    end

    filename = joinpath(DATADIR, "b" * anum[2:end] * ".txt")
    info("Writing " * anum * " to " * filename)

    n = 1
    f = open(filename, "w")
    for l in list
        println(f, n, " ", list[n])
        n += 1
    end
    close(f)
end

doc"""
Make sure that the length of the data section of an OEIS entry does not exceed
260 characters.
"""
function oeis_trimdata(fun, offset::Int)
    len = n = 0
    S = ""
    while true
        st = string(fun(offset + n))
        len += length(st)
        len > 260 && break
        S *= st * ", "
        len += 2
        n += 1
    end
    println(n, " terms")
end

doc"""
Download the sequence with A-number 'anum' from the OEIS to a file in json format.
The file is saved in the 'data' directory.
"""
function oeis_remote(anum)
    if !ismatch(r"^A[0-9]{6}$", anum)
        warn("Not a valid A-number!")
        return
    end

    filename = joinpath(DATADIR, anum * ".json")

    url = URI("http://oeis.org/search?q=id:" * anum * "&fmt=json")
    tries = 3
    r = nothing
    for i = 1:tries
        try
            r = get(url; timeout=.5)
            r.status == 200 && break
            if contains(r.headers["Content-Type"], "text/html")
                display("text/html", r.data)
            end
            r.status == 302 && break # redirection
        catch e
            # warn(e)
        end
        sleep(2)
    end
    if r ≠ nothing && r.status == 200
        open(filename, "w") do f
            write(f, r.data)
        end
    else
        if r == nothing
            warn("Could not download $url, connection timed out.\n")
        else
            warn("Could not download $url\nStatus: $(r.status)")
        end
    end
end

doc"""
Get the sequence with A-number 'anum' from a local copy of the expanded
'stripped.gz' file which can be downloaded from the OEIS. 'bound' is an upper
bound for the number of terms returned. The 'stripped' file is assumed to be in
the 'MathIntSeq/data' directory .
"""
function oeis_local(anum::String, bound::Int)

    if !ismatch(r"^A[0-9]{6}$", anum)
        warn("Not a valid A-number!")
        return []
    end

    oeis_notinstalled() && return []

    A = Array{String}
    data = open(oeis_path())
    for ln in eachline(data)
        if startswith(ln, anum)
            A = split(chop(chomp(ln)), ","; limit=bound + 2)
            break;
        end
    end
    close(data)

    SeqArray([convert(fmpz, parse(BigInt, n)) for n in A[2:min(bound + 1, end)]])
end

doc"""
Search for a sequence in the local OEIS database ('MathIntSeq/data/stripped').
Input the sequence as a comma separated string. If restart = true the search
is redone in the case that no match was found with the first term removed
from the search string. Prints the matches.
"""
function oeis_search(seq::String, restart::Bool)

    oeis_notinstalled() && return []

    found = false
    seq = replace(seq, ' ', "")
    println("Searching for:")
    println(seq)

    data = open(oeis_path())
    for ln in eachline(data)
        index = searchindex(ln, seq)
        index == 0 && continue
        println("Starts at ", index - 10, " ", ln)
        found = true
    end
    close(data)

    if !found && restart
        ind = search(seq, ',')
        if (ind > 0) && (length(seq) > ind)
            seq = seq[ind + 1:end]
            println("Restarting omitting the first term.")
            oeis_search(seq, false)
        end
    end
end

end # module

module OEISUtilsTest
using OEISUtils

function test()

    if oeis_isinstalled()
        info("OEIS data is installed as:")
        info(oeis_path())
    end

    oeis_notinstalled()
end

function demo()
    oeis_writebfile("A000290", n -> n * n, 0, 100)
    println(oeis_trimdata(n -> (-1)^n * n^5, 10))
    oeis_remote("A123456")

    println(oeis_local("A123456", 12))
    println(oeis_local("A015108", 11))

    # A015108 ,1,1,-10,-1231,1636130,23957879562,-3858392581773300,-6835385537899011365535,
    # 133202313157282627679850238250,28553099061411464607955930776882965774,

    seq = "0, 1, 2, 3, 6, 9, 14, 22, 32, 46, 66, 93, 128, 176, 238, 319"
    oeis_search(seq, true)
end

function perf()
end

function main()
    demo()
    test()
    perf()
end

main()

end # module
