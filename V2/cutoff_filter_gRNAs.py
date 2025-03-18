import argparse
import re

def filter_fasta_by_score(input_fasta, output_fasta, score_threshold):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        lines = infile.readlines()
        
        current_header = ""
        current_sequence = ""
        excluded_count = 0  # Counter for excluded entries
        initial_number_entries = 0  # Counter for the total number of entries
        
        for line in lines:
            if line.startswith(">"):  # Header line
                if current_header:
                    # Write the previous gene entry if its score is greater than the threshold
                    match = re.search(r"Sc([\d.]+)", current_header)
                    if match:
                        score = float(match.group(1))
                        if score > score_threshold:
                            outfile.write(current_header + "\n")
                            outfile.write(current_sequence + "\n")
                        else:
                            excluded_count += 1  # Increment the counter for excluded entries
                
                # Start a new gene entry
                current_header = line.strip()
                current_sequence = ""
                initial_number_entries += 1  # Increment the total number of entries
            else:
                current_sequence += line.strip()
        
        # Final gene entry check
        if current_header:
            match = re.search(r"Sc([\d.]+)", current_header)
            if match:
                score = float(match.group(1))
                if score > score_threshold:
                    outfile.write(current_header + "\n")
                    outfile.write(current_sequence + "\n")
                else:
                    excluded_count += 1  # Increment the counter for excluded entries

    return excluded_count, initial_number_entries  # Return the number of excluded and total entries


def main():
    # Argument parser
    parser = argparse.ArgumentParser(description="Filter FASTA file based on score (Sc).")
    parser.add_argument("input_fasta", help="Path to the input FASTA file")
    parser.add_argument("output_fasta", help="Path to the output filtered FASTA file")
    parser.add_argument("score_threshold", type=float, help="Minimum score (Sc) to keep the gene entry")
    
    args = parser.parse_args()
    
    # Call the filter function
    excluded_count, initial_number_entries = filter_fasta_by_score(args.input_fasta, args.output_fasta, args.score_threshold)
    
    print(f"{initial_number_entries-excluded_count} entries were filtered and saved to {args.output_fasta}")
    print(f"{excluded_count} entries with a score < {args.score_threshold} were removed from the total {initial_number_entries}.")

if __name__ == "__main__":
    main()

