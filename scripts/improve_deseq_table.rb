#!/usr/bin/env ruby
### add gene names, biotype, ... to deseq tables

## read in GTF
gtf = {}
File.open(ARGV[0],'r').each do |line|
  s = line.split("\t")
  chr = s[0]
  start = s[3]
  stop = s[4]
  strand = s[6]
  gene_id = s[8].split('gene_id')[1].split(';')[0].gsub('"','').strip
  gene_name = 'NA'
  if line.include?('gene_name')
    gene_name = s[8].split('gene_name')[1].split(';')[0].gsub('"','').strip
  end
  gene_biotype = 'NA'
  if line.include?('gene_biotype')
    gene_biotype = s[8].split('gene_biotype')[1].split(';')[0].gsub('"','').strip
  end

  gtf[gene_id] = [gene_name, gene_biotype, "#{chr}:#{start}-#{stop}", strand]
end
puts "read in #{gtf.size} gtf entries"

#Dir.glob("/mnt/mahlzeitlocal/projects/myotis_rnaseq_weber/rerun01/deseq2/mlu/star/RNA/uniq/*/deseq2*.csv").each do |csv|
Dir.glob("#{ARGV[1]}").each do |csv|
  tmp = File.open(File.dirname(csv) + '/tmp.csv','w')
  csv_file = File.open(csv,'r')
  csv_file.each do |line|
    break if line.include?('gene_id')
    if line.start_with?('""')
      tmp << line.sub('"",','"gene_id","gene_name","gene_biotype","pos","strand",')
    else
      gene_id = line.split(',')[0].gsub('"','')
      tmp << "\"#{gene_id}\"," << gtf[gene_id].join(',') << ',' << line.split(',')[1,8].join(',')
    end
  end
  csv_file.close
  tmp.close

  `mv #{tmp.path} #{csv_file.path}`
  #`ssconvert #{csv_file.path} #{csv_file.path.sub('.csv','.xlsx')}`
  puts csv_file.path
end
