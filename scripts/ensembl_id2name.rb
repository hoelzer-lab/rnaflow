#!/usr/bin/env ruby

gtf = File.open(ARGV[0],'r')
out = File.open(ARGV[1],'w')

gtf.each do |line|
    split = line.split("\t")
    desc = split[8]
    gene_id = desc.split('gene_id')[1].split(';')[0].gsub('"','').chomp.strip
    if line.include?('gene_name')
        gene_name = desc.split('gene_name')[1].split(';')[0].gsub('"','').chomp.strip  
    else
        gene_name = gene_id
    end
    gene_biotype = desc.split('gene_biotype')[1].split(';')[0].gsub('"','').chomp.strip  
    out << "#{gene_id}\t#{gene_name}\t#{gene_biotype}\n"
end

gtf.close
out.close
