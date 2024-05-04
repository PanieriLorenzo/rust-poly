export def yes-or-no [prompt: string] -> bool {
    return ([[name value]; [yes true] [no false]] | input list -d name $prompt).value
}
