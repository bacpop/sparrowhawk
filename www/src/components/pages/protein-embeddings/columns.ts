import { h } from 'vue'
import type { ColumnDef } from '@tanstack/vue-table'
import { ArrowUpDown, ArrowUp, ArrowDown } from "@lucide/vue"
import { Button } from '@/components/ui/button'

/** One row per embedded protein. */
export interface ProteinEmbeddingRow {
    id: string
    length: number
    truncated: boolean
    /** UMAP position, pre-formatted: the pair is read together, never sorted apart. */
    coords: string
}

// Same helper as taxonomic-id/columns.ts; kept local so neither tab's table depends on the
// other's column definitions.
function sortableHeader(label: string, align?: 'right') {
    return ({ column }: { column: { toggleSorting: (asc: boolean) => void; getIsSorted: () => string | false } }) => {
        const sorted = column.getIsSorted()
        const icon = sorted === 'asc' ? ArrowUp : sorted === 'desc' ? ArrowDown : ArrowUpDown
        return h(Button, {
            variant: 'ghost',
            class: align === 'right' ? 'w-full justify-end' : undefined,
            onClick: () => column.toggleSorting(sorted === 'asc'),
        }, () => [label, h(icon, { class: 'ml-2 h-4 w-4' })])
    }
}

export const columns: ColumnDef<ProteinEmbeddingRow>[] = [
    {
        accessorKey: 'id',
        header: sortableHeader('Protein'),
        cell: ({ row }) => h('div', { class: 'font-medium font-mono truncate' }, row.getValue('id')),
    },
    {
        accessorKey: 'length',
        header: sortableHeader('Length', 'right'),
        cell: ({ row }) => h('div', { class: 'text-right' }, `${row.getValue('length')} aa`),
    },
    {
        accessorKey: 'truncated',
        header: sortableHeader('Truncated', 'right'),
        cell: ({ row }) => {
            const t = row.getValue('truncated') as boolean
            return h('div', { class: ['text-right', t ? 'text-yellow-600' : 'text-gray-400'] }, t ? 'yes' : '—')
        },
    },
    {
        // Plain header, not sortable: the cell is formatted text, so a sort would be
        // lexicographic — "9.10, 1.00" ahead of "10.00, 1.00".
        accessorKey: 'coords',
        header: 'UMAP',
        cell: ({ row }) => h('div', { class: 'font-mono truncate' }, row.getValue('coords')),
    },
]
