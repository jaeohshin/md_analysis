from collections import defaultdict

def process_file(filename):
    # Initialize the counts for each pattern (set the order explicitly)
    counts = defaultdict(int)
    
    # Define the correct order of output patterns
    order = [
        "DFGin-ABAminus", "DFGin-BLAminus", "DFGin-BLAplus", "DFGin-BLBminus", 
        "DFGin-BLBplus", "DFGin-BLBtrans", "DFGin-Unassigned", 
        "DFGinter-BABtrans", "DFGinter-Unassigned", "DFGout-BBAminus", 
        "DFGout-Unassigned", "Unassigned-Unassigned"
    ]
    
    # Initialize the dictionary with keys from the order
    for item in order:
        counts[item] = 0
    
    # Read the file and process it line by line
    with open(filename, 'r', encoding='utf-8') as file:
        lines = file.readlines()
        
        # Iterate through every 5 lines (representing one frame)
        for i in range(0, len(lines), 5):
            frame_lines = lines[i:i+5]  # Select the current 5 lines (1 frame)
            
            # Process the 4th column from each of these lines
            frame_info = []
            for line in frame_lines:
                parts = line.split()
                if len(parts) >= 4:
                    fourth_column = parts[3]  # Extract the 4th column (index 3)
                    
                    # Split by underscore and grab the second and third parts
                    split_forth = fourth_column.split('_')
                    
                    # Only consider the second and third parts if they exist
                    if len(split_forth) >= 3:
                        second_part = split_forth[1]
                        third_part = split_forth[2]
                        
                        # Check for 'None' and assign "Unassigned" if found
                        if second_part == 'None' or third_part == 'None':
                            result = "Unassigned-Unassigned"
                        else:
                            result = f"{second_part}-{third_part}"
                    else:
                        result = "Unassigned-Unassigned"
                    
                    # Store the result for this frame
                    frame_info.append(result)
            
            # If we have 5 valid items in the frame_info list, count them
            if len(frame_info) == 5:
                # We count each result occurrence (note we divide by 5 later)
                for result in frame_info:
                    counts[result] += 1
    
    # Divide the counts by 5 to reflect per-frame counts
    for item in order:
        counts[item] = counts[item]  # Each frame has 5 lines, divide the count by 5
    
    # Print the results in the specified order
    for item in order:
        print(f"{item}: {counts[item]}")

# Call the function with your file
filename = "output.txt"  # Replace with the actual file name
process_file(filename)

