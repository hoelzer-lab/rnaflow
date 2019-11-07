#!/usr/bin/env ruby

require 'fileutils'

class RefactorReportingtoolsTable

  def initialize(html_path, anno)

    species = anno.sub('.gene.gtf','')

    case species
      when 'eco' 
        $scan_gene_id_pattern = 'ER[0-9]+_[0-9]+'
        $ensembl_url = 'https://bacteria.ensembl.org/Escherichia_coli_k_12/Gene/Summary?g='
      when 'hsa'
        $scan_gene_id_pattern = 'ENSG[0-9]+'
        $ensembl_url = 'https://ensembl.org/Homo_sapiens/Gene/Summary?g='
      else
        $scan_gene_id_pattern = false
        $ensembl_url = false
    end

		$id2name = {}
    $id2biotype = {}
    $id2pos = {}
    gtf = File.open(anno,'r')
    gtf.each do |line|
      s = line.split("\t")
      if line.include?('gene_name')
        gene_name = s[8].split('gene_name')[1].split(';')[0].gsub('"','').strip
      else
        gene_name = 'NA'
      end
      gene_id = s[8].split('gene_id')[1].split(';')[0].gsub('"','').strip
      gene_biotype = s[8].split('gene_biotype')[1].split(';')[0].gsub('"','').strip
      chr = s[0]
      start = s[3]
      stop = s[4]
      strand = s[6]
      $id2name[gene_id] = gene_name
      $id2biotype[gene_id] = gene_biotype
      $id2pos[gene_id] = [chr, start, stop, strand]
    end
    gtf.close
    puts "read in #{$id2name.keys.size} genes."

    refactor_deseq_html_table(html_path)
  end
  
  
  def refactor_deseq_html_table(html_path)
    
    html_file = File.open(html_path,'r')    
    html_file_refac = File.open(html_path.sub('.html','_table.html'),'w')

    pvalue = File.basename(html_path, '.html').split('_').reverse[0]

    html_file.each do |line|
      if line.start_with?('<div class="container"')
        tmp_array = line.split('</tr>')

        html_file_refac << refac_table_header(tmp_array[0], 1) << '</tr>' # first header
        html_file_refac << refac_table_header(tmp_array[1], 2) << '</tr>' # second header
        tmp_array[2..tmp_array.length-3].each do |row|
          row_splitted = row.sub('<a href','<a target="_blank" href').split('</td>')
          gene_id = row_splitted[0].scan(/#{$scan_gene_id_pattern}/)[0]
          gene_id = gene_id.gsub('"','')
          gene_name = $id2name[gene_id]
          gene_biotype = $id2biotype[gene_id]
          #puts gene_id          
          pos_part = "<td class=\"\">#{$id2pos[gene_id][0]}:#{$id2pos[gene_id][1]}-#{$id2pos[gene_id][2]} (#{$id2pos[gene_id][3]})"
          if gene_id.include?('ER')
	          new_row = [row_splitted[0].sub('<td class="">',"<td class=\"\"><a target=\"_blank\" href=\"#{$ensembl_url}#{gene_id};\">") + '</a>', "<td class=\"\">#{gene_name}", "<td class=\"\">#{gene_biotype}", pos_part, row_splitted[1].sub('href=','target="_blank" href=').sub('<td class="">','<td class=""><div style="width: 200px">') + '</div>', row_splitted[2], row_splitted[3], row_splitted[4]].join('</td>') << '</td>'
	        else
            new_row = [row_splitted[0], "<td class=\"\">#{gene_name}", "<td class=\"\">#{gene_biotype}", pos_part, row_splitted[1].sub('href=','target="_blank" href=').sub('<td class="">','<td class=""><div style="width: 200px">') + '</div>', row_splitted[2], row_splitted[3], row_splitted[4]].join('</td>') << '</td>'            
          end          
          html_file_refac << new_row << '</tr>'
          
        end
        table_footer = tmp_array[tmp_array.length-2]
        html_file_refac << refac_table_header(table_footer, 3) << '</tr>' # table footer
        html_file_refac << tmp_array[tmp_array.length-1] # footer
      else
        html_file_refac << line
      end
    end
    html_file.close; html_file_refac.close
  end


  def refac_table_header(string, type)    
    split = string.split('</th>')    
    case type
      when 1
        [split[0], "<th class=\" sort-string-robust top-header-row no-print\">Name", "<th class=\" sort-string-robust top-header-row no-print\">Type", "<th class=\" sort-string-robust top-header-row no-print\">Position", split[1], split[2], split[3], split[4]].join('</th>') << '</th>'
      when 2
        [split[0], "<th class=\" sort-string-robust bottom-header-row\">Name", "<th class=\" sort-string-robust bottom-header-row\">Type", "<th class=\" sort-string-robust bottom-header-row\">Position", split[1], split[2], split[3], split[4]].join('</th>') << '</th>'
      else
        [split[0], "<th class=\"sort-string-robust bottom-header-row\">Name", "<th class=\"sort-string-robust bottom-header-row\">Type", "<th class=\"sort-string-robust bottom-header-row\">Position", split[1], split[2], split[3], split[4]].join('</th>') << '</th>'
    end
  end

end

RefactorReportingtoolsTable.new(ARGV[0], ARGV[1])
