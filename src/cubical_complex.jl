function parse_im_data(file_path)
    fid = open(file_path, "r")

    # Skip the 1024-byte header
    seek(fid, 1024)  # Move the file pointer to the first byte after the header

    hex_array = read(fid)
    integer_array = [Int(hex) for hex in hex_array]
    data3d = reshape(integer_array, (128, 128, 128))  

    return(data3d)
end

file_path = "resources/snakesIm/1.im"
parse_im_data(file_path)
