# C6C-3 provider contact verification (public pages only; nothing sent)

Status vocabulary: VERIFIED (same address on >=2 independent official pages), LISTED (one official page), UNRESOLVED.

| Source | Role | Route | Evidence | Status |
|---|---|---|---|---|
| GPPVA_FSU2 | Data-sheet submitter: Helena Pais (Univ. Coimbra, CFisUC) | `hpais@uc.pt` | Data sheet `Crust/FSU2/eos.pdf` ("submitted by Helena Pais", "e-mail contact hpais@uc.pt") and CFisUC profile https://cfisuc.fis.uc.pt/people.php?oid=5043442 | **VERIFIED** |
| SPG_M1 | Data-sheet submitter: Luigi Scurto | `lscurto@student.uc.pt` | Data sheet `Crust/M1/eos.pdf` only. The CFisUC profile (https://cfisuc.fis.uc.pt/people.php?oid=12300255) describes him as "Former PhD Student", so a student address may be stale. | LISTED (possibly stale) |
| SPG_M1 | same person, institutional page | `luigiscurto.pr@gmail.com` | Listed on the CFisUC profile above. Not on the data sheet; personal-domain address published by the institution. | LISTED |
| SPG_M1 | Senior co-author / group lead of the matching paper: Helena Pais (also PRD 109, 103015 author) | `hpais@uc.pt` | Same verified route as above; a reasonable copy recipient so the question reaches someone with access to the export pipeline. Whether she ran the SPG export is **not established**. | VERIFIED route, role unconfirmed |
| Both | CompOSE core team | address shown on https://compose.obspm.fr/contacts as "you can reach this team of core members under the following e-mail address", rendered through a client-side email obfuscation. | Page fetched (HTTP 200). The address was deliberately not decoded or transcribed; a search summary reported a `develop.*@obspm.fr` address, which is unverified here. Use the address as displayed in a browser. | LISTED (not transcribed) |

Not found/unresolved: the person who actually ran the CompOSE export for each product (the data sheets state only the submitter). Gulminelli (co-author) was not looked up; not required.
Not used: personal/private contacts, scraped directories, phone numbers. `www.uc.pt` returned HTTP 403 to a plain fetch and was not pursued.

Suggested routing (for the project owner to decide): GPPVA_FSU2 -> Pais; SPG_M1 -> Scurto (CFisUC-listed address), cc Pais; optionally cc CompOSE core team for metadata/data-sheet clarification.
