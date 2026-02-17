import os
import re
from dataclasses import dataclass
from rdkit import Chem
from rdkit import RDLogger

RDLogger.EnableLog("rdApp.error")


@dataclass
class Smarts:
    smarts: str = ""
    rawsmarts: str = ""
    name: str = ""
    groupname: str = ""
    search: Chem.Mol | None = None

    def __str__(self):
        return f"{self.smarts} {self.name}".strip()


class SmartsFile:
    def __init__(self):
        self.smartses: list[Smarts] = []
        self.defines: dict[str, str] = {}
        self.failed_smarts: list[dict] = []
        self.groupnames: list[str] = []
        self.rawtxt: str = ""
        self.name: str = ""

    def size(self) -> int:
        return len(self.smartses)

    def get_smarts_object(self, i: int) -> Smarts:
        return self.smartses[i]

    def parse_file(self, ftxt_or_path: str, strict: bool = False, groupname: str = "") -> bool:
        if os.path.isfile(ftxt_or_path):
            self.name = os.path.basename(ftxt_or_path)
            with open(ftxt_or_path, "r", encoding="utf-8", errors="replace") as f:
                ftxt = f.read()
        else:
            ftxt = ftxt_or_path

        p_define = re.compile(r"^\s*define\s+\$(\w+)\s+(.+?)\s*$", re.IGNORECASE)
        p_smarts_name = re.compile(r"^\s*(\S+)(?:\s+(.+))?\s*$")

        for line_num, raw_line in enumerate(ftxt.splitlines(), 1):
            self.rawtxt += raw_line + "\n"
            line = self._strip_inline_comment(raw_line).strip()
            if not line:
                continue

            m = p_define.match(line)
            if m:
                tag = m.group(1)
                val = m.group(2).strip()
                self.defines[tag] = val
                continue

            m2 = p_smarts_name.match(line)
            if not m2:
                continue

            rawsmarts = (m2.group(1) or "").strip()
            name = (m2.group(2) or "").strip()

            try:
                expanded = self._expand_all_macros(rawsmarts, strict=strict)

                self._pre_validate_smarts(expanded)

                fixed = self._normalize_for_rdkit(expanded)

                self._post_validate_smarts(fixed)

                q = Chem.MolFromSmarts(fixed)
                if q is None:
                    fixed = self._attempt_recovery(fixed)
                    q = Chem.MolFromSmarts(fixed)
                    if q is None:
                        raise ValueError(f"RDKit could not parse: {fixed}")

                if groupname and groupname not in self.groupnames:
                    self.groupnames.append(groupname)

                self.smartses.append(Smarts(
                    smarts=fixed,
                    rawsmarts=rawsmarts,
                    name=name,
                    groupname=groupname,
                    search=q
                ))

            except Exception as e:
                expanded_str = expanded if 'expanded' in dir() else None
                fixed_str = fixed if 'fixed' in dir() else None
                
                if strict:
                    raise Exception(
                        f"Problem parsing SMARTS (line {line_num}):\n"
                        f"  raw: {rawsmarts}\n"
                        f"  expanded: {expanded_str}\n"
                        f"  fixed: {fixed_str}\n"
                        f"  name: {name}\n"
                        f"  error: {e}"
                    ) from e

                self.failed_smarts.append({
                    "line": line_num,
                    "raw": rawsmarts,
                    "name": name,
                    "expanded": expanded_str,
                    "fixed": fixed_str,
                    "error": str(e),
                })

        for i, s in enumerate(self.smartses, 1):
            if not s.name:
                s.name = str(i)

        return True
    _tag_re = re.compile(r"\$([A-Za-z0-9_]+)")
    _bare_macro_re = re.compile(r"(?<!\$)\b([A-Z][A-Za-z0-9_]{2,})\b")

    def _expand_all_macros(
        self,
        smarts: str,
        strict: bool,
        max_rounds: int = 200,
    ) -> str:
        seen_states = set()
        cur = smarts

        for _ in range(max_rounds):
            if cur in seen_states:
                raise ValueError("Circular macro expansion detected.")
            seen_states.add(cur)

            m = self._tag_re.search(cur)
            if not m:
                return cur 

            new_cur = self._smart_expand_macros(cur)

            if new_cur == cur:
                m = self._tag_re.search(cur)
                if m and m.group(1) not in self.defines:
                    if strict:
                        raise ValueError(f"Missing define: {m.group(1)}")
                    break
                cur = self._tag_re.sub(lambda mm: self._replace_tag(mm, strict), cur)
            else:
                cur = new_cur

        return cur

    def _replace_tag(self, match: re.Match, strict: bool) -> str:
        tag = match.group(1)
        if tag not in self.defines:
            if strict:
                raise ValueError(f"Missing define: {tag}")
            return match.group(0)  

        return self.defines[tag]

    def _smart_expand_macros(self, smarts: str) -> str:
        result = []
        i = 0
        while i < len(smarts):
            m = self._tag_re.match(smarts, i)
            if m:
                tag = m.group(1)
                if tag in self.defines:
                    replacement = self.defines[tag]
                    prefix = smarts[:i]
                    in_bracket = prefix.count('[') - prefix.count(']') > 0
                    
                    is_molecule = self._is_molecule_pattern(replacement)
                    is_atom_list = replacement.startswith('[') and replacement.endswith(']') and not is_molecule
                    
                    if in_bracket:
                        if is_atom_list:
                            replacement = replacement[1:-1]
                        elif is_molecule:
                            replacement = f'$({replacement})'
                    
                    result.append(replacement)
                    i = m.end()
                else:
                    result.append(smarts[i])
                    i += 1
            else:
                # Check for bare macro name without $ prefix (e.g. !B8EXC → !$B8EXC)
                bm = self._bare_macro_re.match(smarts, i)
                if bm and bm.group(1) in self.defines:
                    tag = bm.group(1)
                    replacement = self.defines[tag]
                    prefix = smarts[:i]
                    in_bracket = prefix.count('[') - prefix.count(']') > 0

                    is_molecule = self._is_molecule_pattern(replacement)
                    is_atom_list = replacement.startswith('[') and replacement.endswith(']') and not is_molecule

                    if in_bracket:
                        if is_atom_list:
                            replacement = replacement[1:-1]
                        elif is_molecule:
                            replacement = f'$({replacement})'

                    result.append(replacement)
                    i = bm.end()
                else:
                    result.append(smarts[i])
                    i += 1
        return ''.join(result)

    def _is_molecule_pattern(self, pattern: str) -> bool:
        if not pattern:
            return False
        
        if pattern.startswith('[') and pattern.endswith(']'):
            # Verify these are matching brackets (i.e. not [A][B] or [A](B))
            is_pair = True
            depth = 0
            for i, ch in enumerate(pattern[:-1]):
                if ch == '[':
                    depth += 1
                elif ch == ']':
                    depth -= 1
                
                if depth == 0 and i > 0:
                    is_pair = False
                    break
            
            if is_pair:
                inner = pattern[1:-1]
                bracket_depth = 0
                for ch in inner:
                    if ch == '[':
                        bracket_depth += 1
                    elif ch == ']':
                        bracket_depth -= 1
                    elif bracket_depth == 0:
                        if ch in '=#~:-' or ch.isdigit():
                            return True
                        if ch == '(':
                            return True
                return False

        bracket_depth = 0
        for ch in pattern:
            if ch == '[':
                bracket_depth += 1
            elif ch == ']':
                bracket_depth -= 1
            elif bracket_depth == 0:
                if ch in '=#~:-()' or ch.isdigit():
                    return True
        
        atom_count = sum(1 for ch in pattern if ch.isupper() and ch in 'CNOSPFIBKH')
        return atom_count > 1 or any(ch in pattern for ch in '()=#~:-')

    def _strip_inline_comment(self, line: str) -> str:
        in_brackets = 0
        for i, ch in enumerate(line):
            if ch == "[":
                in_brackets += 1
            elif ch == "]" and in_brackets > 0:
                in_brackets -= 1
            elif ch == "#" and in_brackets == 0:
                if i == 0:
                    return ""  
                prev_char = line[i-1]
                if prev_char.isspace():
                    return line[:i].rstrip()
        return line

    def _strip_single_letter_parens(self, smarts: str) -> str:
        """Strip single-letter parens like (H) or (C), but NOT inside $() recursive SMARTS."""
        result = []
        i = 0
        dollar_depth = 0  # track depth inside $(...) 
        while i < len(smarts):
            if smarts[i:i+2] == '$(':
                dollar_depth += 1
                result.append('$(')
                i += 2
                continue
            if dollar_depth > 0:
                if smarts[i] == '(':
                    dollar_depth += 1
                elif smarts[i] == ')':
                    dollar_depth -= 1
                result.append(smarts[i])
                i += 1
                continue
            # Outside $(): apply the (X) → X rule
            if (smarts[i] == '(' and i + 2 < len(smarts) 
                    and smarts[i+2] == ')' and smarts[i+1].isalpha()):
                result.append(smarts[i+1])
                i += 3
            else:
                result.append(smarts[i])
                i += 1
        return ''.join(result)

    def _fix_bare_ch(self, smarts: str) -> str:
        """Convert bare 'CH' outside brackets to '[CH]' in SMARTS notation.
        
        In SMARTS, bare 'CH' means C followed by H (two atoms). But in many .sma
        files, 'CH' means carbon with one implicit hydrogen, which is '[CH]'.
        Only convert when CH is followed by a bond, bracket, or ring-start character.
        """
        result = []
        i = 0
        bracket_depth = 0
        dollar_depth = 0
        while i < len(smarts):
            if smarts[i] == '[':
                bracket_depth += 1
            elif smarts[i] == ']':
                bracket_depth -= 1
            if smarts[i:i+2] == '$(':
                dollar_depth += 1
            if dollar_depth > 0 and smarts[i] == ')':
                dollar_depth -= 1
                
            # Only convert CH outside brackets and outside $()
            if (bracket_depth == 0 and dollar_depth == 0
                    and i + 1 < len(smarts)
                    and smarts[i] == 'C' and smarts[i+1] == 'H'
                    and (i + 2 >= len(smarts) or smarts[i+2] in '=~-:#@!()[]0123456789')):
                # Check it's not already inside brackets or following another letter
                if i == 0 or smarts[i-1] in '=~-:#@!()[]0123456789,;>':
                    result.append('[CH]')
                    i += 2
                    continue
            result.append(smarts[i])
            i += 1
        return ''.join(result)

    def _fix_bare_recursive_smarts(self, smarts: str) -> str:
        s = smarts
        
        if s.startswith('$(') and not s.startswith('[$('):
            depth = 0
            for i, ch in enumerate(s):
                if ch == '(':
                    depth += 1
                elif ch == ')':
                    depth -= 1
                    if depth == 0:
                        if i == len(s) - 1:
                            s = '[' + s + ']'
                        break
        
        return s

    def _fix_and_symbol(self, smarts: str) -> str:
        result = []
        bracket_depth = 0
        
        for ch in smarts:
            if ch == '[':
                bracket_depth += 1
                result.append(ch)
            elif ch == ']':
                bracket_depth -= 1
                result.append(ch)
            elif ch == '&' and bracket_depth > 0:
                result.append(';')
            else:
                result.append(ch)
        
        return ''.join(result)

    def _fix_nested_brackets(self, smarts: str) -> str:
        def fix_nested(match):
            outer_atom = match.group(1)
            inner_list = match.group(2)
            if '$' not in outer_atom:
                return f'[{outer_atom}]~[{inner_list}]'
            return match.group(0)

        s = re.sub(r'\[([A-Za-z][A-Za-z0-9;,+\-]*)\[([A-Za-z][A-Za-z0-9,;+\-]*)\]\]', fix_nested, smarts)
        
        return s

    def _fix_molecule_in_brackets(self, smarts: str) -> str:
        return smarts

    def _contains_molecule_structure(self, pattern: str) -> bool:
        if not pattern:
            return False
            
        bracket_depth = 0
        paren_depth = 0
        
        for i, ch in enumerate(pattern):
            if ch == '[':
                bracket_depth += 1
            elif ch == ']':
                bracket_depth -= 1
            elif ch == '(':
                paren_depth += 1
            elif ch == ')':
                paren_depth -= 1
            elif bracket_depth == 0:
                if ch in '=#~:':
                    return True
                if ch.isdigit() and i > 0 and pattern[i-1] not in '+-':
                    return True
                if ch == '-' and i > 0 and pattern[i-1] not in '[-;' and (i+1 < len(pattern) and pattern[i+1] not in '0123456789]'):
                    return True
        
        return False

    def _split_on_comma_outside_brackets(self, s: str) -> list:
        parts = []
        current = []
        bracket_depth = 0
        paren_depth = 0
        
        for ch in s:
            if ch == '[':
                bracket_depth += 1
            elif ch == ']':
                bracket_depth -= 1
            elif ch == '(':
                paren_depth += 1
            elif ch == ')':
                paren_depth -= 1
            
            if ch == ',' and bracket_depth == 0 and paren_depth == 0:
                parts.append(''.join(current))
                current = []
            else:
                current.append(ch)
        
        if current:
            parts.append(''.join(current))
        
        return parts

    def _fix_mixed_bracket_content(self, smarts: str) -> str:
        return smarts
    def _strip_smastar_labels(self, smarts: str) -> str:
        return re.sub(r'_(\d+)\b', '', smarts)

    def _normalize_for_rdkit(self, smarts: str) -> str:
        s = smarts
        s = self._strip_smastar_labels(s)
        s = self._strip_outer_parens_if_wrapping_whole(s)
        s = self._fix_dash_star_branches(s)
        s = self._fix_dollar_star_atom_shorthand(s)

        s = self._fix_and_symbol(s)
        
        s = re.sub(r'\(\*\)', '*', s)
        # Strip single-letter parens like (H) or (C), but NOT inside $() recursive SMARTS
        s = self._strip_single_letter_parens(s)
        
        s = self._fix_nested_brackets(s)

        s = self._fix_molecule_in_brackets(s)
        
        s = re.sub(r'([=#~\-:]),(?=[\]\)]|$)', r'\1', s)
        s = re.sub(r'([=#~\-:]),$', r'\1', s)
        s = re.sub(r'=,\)', '=)', s) 
        
        s = re.sub(r'\$\(([+\-]\d)', r'$([\1', s)

        s = re.sub(r'\$\(\$\(', '$(', s)

        s = self._fix_bare_recursive_smarts(s)

        # Convert bare CH outside brackets to [CH]
        s = self._fix_bare_ch(s)

        s = s.strip()
        
        s = re.sub(r'\s+(?=[)\]])', '', s)  
        s = re.sub(r'(?<=[([])\s+', '', s)  

        s = re.sub(r"S\(\(=O\)=O\)", r"S(=O)(=O)", s)

        s = re.sub(r"C\(\(F\)F\)F", r"C(F)(F)F", s)
        
        s = re.sub(r"\(\(([A-Za-z])\)([A-Za-z])\)", r"(\1)(\2)", s)
        
        s = re.sub(r"\(\(\s*(\[[^\]]+\])\s*\)\s*(\[[^\]]+\])\s*\)", r"(\1\2)", s)

        s = re.sub(r'\[\]', '[*]', s)
        
        s = re.sub(r'\[,([^\]]+)\]', r'[\1]', s)
        s = re.sub(r'\[([^\]]+),\]', r'[\1]', s)

        s = re.sub(r',{2,}', ',', s)

        s = re.sub(r'\[;([^\]]+)\]', r'[\1]', s)
        s = re.sub(r'\[([^\]]+);\]', r'[\1]', s)

        s = re.sub(r'(?<![!\$])={2,}(?!=)', '=', s)
        
        s = re.sub(r'(?<![!\$])-{2,}(?!>)', '-', s)
        
        s = re.sub(r'#{2,}', '#', s)

        s = re.sub(r'%(\d)(?!\d)', r'\1', s)

        s = re.sub(r'\$\(\s*\)', '', s)
        
        s = self._balance_recursive_smarts(s)

        s = re.sub(r'!!', '', s)
        
        s = re.sub(r'\[([^\]]*?)!;([^\]]*?)\]', r'[\1!\2]', s)
        s = re.sub(r'\[([^\]]*?);!([^\]]*?)\]', r'[\1!\2]', s)

        s = re.sub(r'\[([^\]]*?)\+{3,}\]', lambda m: f'[{m.group(1)}+{m.group(0).count("+")}]', s)
        s = re.sub(r'\[([^\]]*?)-{3,}\]', lambda m: f'[{m.group(1)}-{m.group(0).count("-")}]', s)

        s = re.sub(r'\[H-\d+\]', '[H0]', s)

        s = re.sub(r'\.{2,}', '.', s)
        
        s = re.sub(r'^\.+', '', s)
        s = re.sub(r'\.+$', '', s)

        s = re.sub(r'\[0([A-Za-z])', r'[\1', s)

        s = re.sub(r'([=\-#~@:]),$', r'\1', s)
        s = re.sub(r',([=\-#~@:])', r'\1', s)
        s = self._balance_parentheses(s)
        
        s = self._balance_brackets(s)

        return s

    def _balance_recursive_smarts(self, smarts: str) -> str:
        result = []
        i = 0
        while i < len(smarts):
            if smarts[i:i+2] == '$(':
                # Find the matching closing paren
                result.append('$(')
                i += 2
                depth = 1
                start = i
                while i < len(smarts) and depth > 0:
                    if smarts[i] == '(':
                        depth += 1
                    elif smarts[i] == ')':
                        depth -= 1
                    if depth > 0:
                        result.append(smarts[i])
                    i += 1
                if depth > 0:
                    result.append(')')
                else:
                    result.append(')')
            else:
                result.append(smarts[i])
                i += 1
        return ''.join(result)

    def _balance_parentheses(self, smarts: str) -> str:
        depth = 0
        in_bracket = 0
        
        for ch in smarts:
            if ch == '[':
                in_bracket += 1
            elif ch == ']':
                in_bracket = max(0, in_bracket - 1)
            elif ch == '(' and in_bracket == 0:
                depth += 1
            elif ch == ')' and in_bracket == 0:
                depth -= 1
        
        if depth > 0:
            smarts += ')' * depth
        
        return smarts

    def _balance_brackets(self, smarts: str) -> str:
        depth = 0
        
        for ch in smarts:
            if ch == '[':
                depth += 1
            elif ch == ']':
                depth -= 1
        
        if depth > 0:
            smarts += ']' * depth
        
        return smarts
    def _strip_outer_parens_if_wrapping_whole(self, smarts: str) -> str:
        """
        If the entire SMARTS is wrapped in a single redundant pair of parentheses,
        remove them: '(A.B)' -> 'A.B'
        This is safe only if the outer '(' matches the final ')'.
        """
        s = smarts.strip()
        if not (s.startswith('(') and s.endswith(')')):
            return s

        depth = 0
        in_bracket = 0
        for i, ch in enumerate(s):
            if ch == '[':
                in_bracket += 1
            elif ch == ']':
                in_bracket = max(0, in_bracket - 1)
            elif in_bracket == 0:
                if ch == '(':
                    depth += 1
                elif ch == ')':
                    depth -= 1
                    if depth == 0 and i != len(s) - 1:
                        return s

        return s[1:-1].strip()

    def _fix_dash_star_branches(self, smarts: str) -> str:
        s = smarts

        s = re.sub(r'\(\s*-\s*\*\s*\)', r'([*])', s)

        s = re.sub(r'-(\*)', r'-[*]', s)

        return s

    def _pre_validate_smarts(self, smarts: str) -> None:
        if not smarts:
            raise ValueError("Empty SMARTS string")
        
        invalid_chars = re.findall(r'[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]', smarts)
        if invalid_chars:
            raise ValueError(f"Invalid control characters in SMARTS: {invalid_chars}")

    def _post_validate_smarts(self, smarts: str) -> None:
        depth = 0
        in_bracket = 0
        for i, ch in enumerate(smarts):
            if ch == '[':
                in_bracket += 1
            elif ch == ']':
                in_bracket = max(0, in_bracket - 1)
            elif ch == '(' and in_bracket == 0:
                depth += 1
            elif ch == ')' and in_bracket == 0:
                depth -= 1
                if depth < 0:
                    raise ValueError(f"Unmatched closing parenthesis at position {i}")
        
        if depth != 0:
            raise ValueError(f"Unbalanced parentheses: {depth} unclosed")
        
        bracket_depth = 0
        for i, ch in enumerate(smarts):
            if ch == '[':
                bracket_depth += 1
            elif ch == ']':
                bracket_depth -= 1
                if bracket_depth < 0:
                    raise ValueError(f"Unmatched closing bracket at position {i}")
        
        if bracket_depth != 0:
            raise ValueError(f"Unbalanced brackets: {bracket_depth} unclosed")

    def _attempt_recovery(self, smarts: str) -> str:
        s = smarts
        
        if '..' in s:
            s = re.sub(r'\.{2,}', '.', s)
        
        if s.count('$(') > 5: 
            s = re.sub(r'\$\(\$\(([^)]+)\)\)', r'$(\1)', s)
        
        def simplify_negations(match):
            content = match.group(1)
            negations = re.findall(r'![A-Za-z#][0-9]*', content)
            if len(negations) > 10:
                return '[' + ';'.join(negations[:5]) + ']'
            return match.group(0)
        
        s = re.sub(r'\[([^\]]+)\]', simplify_negations, s)
        
        s = re.sub(r'\[a([,;])', r'[c,n,o,s\1', s)
        
        def fix_ring_number(match):
            num = int(match.group(1))
            if num > 99:
                return str(num % 99 + 1)
            return match.group(0)
        
        s = re.sub(r'%(\d{3,})', fix_ring_number, s)
        
        s = re.sub(r'(\d)([=\-#:~])$', r'\1', s)

        s = re.sub(r'!@;=', '!@=', s)
        s = re.sub(r'!@;-', '!@-', s)
        s = re.sub(r'!@;#', '!@#', s)
        s = re.sub(r'!@;~', '!@~', s)

        s = re.sub(r'^[=\-#~:]', '', s)
        s = re.sub(r'[=\-#~:]$', '', s)
        
        s = re.sub(r'\(\)', '', s)
        
        s = re.sub(r'\[#0\]', '[*]', s)

        s = re.sub(r'\$\$\(', '$(', s)
        
        s = re.sub(r'\[\:(\d+)\]', r'[*:\1]', s)
        
        return s

    def print_failed_smarts(self):
        if not self.failed_smarts:
            print("No failures.")
            return
        print("\nFAILED SMARTS:\n")
        for f in self.failed_smarts:
            print(f"Line {f['line']}: {f['raw']}")
            if f.get("name"):
                print(f"  Name: {f['name']}")
            if f.get("expanded") and f["expanded"] != f["raw"]:
                print(f"  Expanded: {f['expanded']}")
            if f.get("fixed") and f["fixed"] != f.get("expanded"):
                print(f"  Fixed: {f['fixed']}")
            print(f"  Error: {f['error']}\n")

    def print_smarts(self, file=None):
        """Print all SMARTS patterns with their names."""
        import sys
        out = file if file else sys.stdout
        for s in self.smartses:
            print(f"{s.smarts} {s.name}", file=out)

    def get_validated_smarts(self) -> list[tuple[str, str, bool]]:
        """
        Get all SMARTS patterns with RDKit validation status.
        Returns list of (smarts, name, is_valid) tuples.
        """
        results = []
        for s in self.smartses:
            is_valid = s.search is not None
            if is_valid:
                try:
                    q = Chem.MolFromSmarts(s.smarts)
                    is_valid = q is not None
                except Exception:
                    is_valid = False
            results.append((s.smarts, s.name, is_valid))
        return results

    def _fix_dollar_star_atom_shorthand(self, smarts: str) -> str:
        s = smarts
        s = re.sub(r'\[H\s*(#\d+)\s*\]', r'[H,\1]', s)
        s = re.sub(r'\[H\s*(#\d+)\s*\]', r'[H,\1]', s)

        s = re.sub(r'\[H\s*([A-Z][a-z]?)\s*\]', r'[H,\1]', s)

        s = re.sub(r'\$\(\*\s*([A-Za-z])\s*\)', r'$(*[\1])', s)

        s = re.sub(r'\$\(\*\[C\]\)', r'$(*[#6])', s)

        # NOTE: _rewrite_bracket_orlist was removed because it incorrectly
        # collapsed recursive SMARTS like [$(*[H]),$(*[#6])] into [H,#6].
        # These are semantically different: $(*[H]) means "any atom bonded to H"
        # which is NOT the same as just [H].

        return s


    def export_to_file(self, output_path: str, include_validation: bool = True) -> int:
        valid_count = 0
        with open(output_path, 'w', encoding='utf-8') as f:
            for s in self.smartses:
                is_valid = s.search is not None
                if is_valid:
                    try:
                        q = Chem.MolFromSmarts(s.smarts)
                        is_valid = q is not None
                    except Exception:
                        is_valid = False
                
                name = s.name.replace(' ', '_')
                if include_validation and not is_valid:
                    # Skip invalid patterns to avoid comments in output
                    pass
                else:
                    f.write(f"{s.smarts} {name}\n")
                    valid_count += 1
        return valid_count


def main():
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(
        description="Parse .sma files and output RDKit-validated SMARTS patterns with names"
    )
    parser.add_argument(
        "sma_file",
        help="Path to the .sma file to parse"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output .txt file path (default: <input_name>_validated.txt)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Show parsing statistics and any failures"
    )
    parser.add_argument(
        "-s", "--strict",
        action="store_true",
        help="Strict mode - fail on any parsing error"
    )
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="Skip RDKit validation (include all patterns)"
    )
    parser.add_argument(
        "--stdout",
        action="store_true",
        help="Print to stdout instead of file"
    )
    
    args = parser.parse_args()
    
    if not os.path.isfile(args.sma_file):
        print(f"Error: File not found: {args.sma_file}", file=sys.stderr)
        sys.exit(1)
    
    sf = SmartsFile()
    
    try:
        sf.parse_file(args.sma_file, strict=args.strict)
    except Exception as e:
        print(f"Error parsing file: {e}", file=sys.stderr)
        sys.exit(1)
    
    if args.output:
        output_path = args.output
    else:
        base_name = os.path.splitext(os.path.basename(args.sma_file))[0]
        output_path = f"{base_name}_validated.txt"
    
    if args.stdout:
        validated = sf.get_validated_smarts()
        valid_count = 0
        for smarts, name, is_valid in validated:
            encoded_name = name.replace(' ', '_')
            if args.no_validate or is_valid:
                print(f"{smarts} {encoded_name}")
                valid_count += 1
            # Skip outputting invalid patterns if verbose, to avoid comments
    else:
        valid_count = sf.export_to_file(output_path, include_validation=not args.no_validate)
        print(f"Exported {valid_count} valid SMARTS patterns to: {output_path}")
    
    if args.verbose:
        print(f"\nStatistics:", file=sys.stderr)
        print(f"  Total parsed: {sf.size()}", file=sys.stderr)
        print(f"  Failed to parse: {len(sf.failed_smarts)}", file=sys.stderr)
        if sf.failed_smarts:
            print(f"\nFailed patterns:", file=sys.stderr)
            sf.print_failed_smarts()


if __name__ == "__main__":
    main()
