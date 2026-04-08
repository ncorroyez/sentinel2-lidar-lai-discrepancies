def vigenere_decrypt(ciphertext, key):
    decrypted = ""
    key = key.upper().replace(" ", "")
    ciphertext = ciphertext.upper().replace(" ", "")
    
    for i in range(len(ciphertext)):
        # Standard Vigenere formula: (C - K) % 26
        char_val = (ord(ciphertext[i]) - ord(key[i % len(key)])) % 26
        decrypted += chr(char_val + ord('A'))
        
    return decrypted

# Full ciphertext and the key that fits the 28 characters
cipher = "NSDWZMNNYLELRFSGIGRFXLJFMUBW"
key = "CESTLAFAUTEDEMONA" # This key repeats to decrypt the message

print(f"Decrypted message: {vigenere_decrypt(cipher, key)}")